//! Progress bar and line counting utilities.
//!
//! Provides real-time progress display with throughput and ETA estimation,
//! plus efficient file line counting for progress tracking.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};

use anyhow::Result;

/// Estimate the number of lines in a file by sampling the first few lines
/// and extrapolating from total file size.
///
/// This avoids a full file scan before processing begins, saving I/O time
/// for large sequence files. The estimate uses the average line length of
/// the first `sample_lines` non-empty lines.
///
/// Returns (estimated_line_count, file_size_bytes).
pub fn estimate_total_lines(path: &str, sample_lines: usize) -> Result<(u64, u64)> {
    let file_size = std::fs::metadata(path)?.len();
    if file_size == 0 {
        return Ok((0, 0));
    }

    let mut file = File::open(path)?;
    let mut buf = [0u8; 256 * 1024];
    let mut total_bytes = 0u64;
    let mut sampled_lines = 0usize;

    'outer: loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        for &b in &buf[..n] {
            total_bytes += 1;
            if b == b'\n' {
                sampled_lines += 1;
                if sampled_lines >= sample_lines {
                    break 'outer;
                }
            }
        }
    }

    // If we sampled enough lines, extrapolate
    if sampled_lines > 0 && total_bytes > 0 {
        let avg_bytes_per_line = total_bytes as f64 / sampled_lines as f64;
        let estimated = (file_size as f64 / avg_bytes_per_line).ceil() as u64;
        file.seek(SeekFrom::Start(0))?;
        Ok((estimated.max(1), file_size))
    } else {
        // File is too short to sample — count exactly
        file.seek(SeekFrom::Start(0))?;
        let exact = count_lines_in_buf(&mut file)?;
        file.seek(SeekFrom::Start(0))?;
        Ok((exact, file_size))
    }
}

/// Count newlines in a file by scanning raw bytes (256KB buffer).
///
/// Fast line counting that avoids UTF-8 validation overhead of
/// line-by-line reading. Handles files with or without trailing newline.
/// Seeks back to the beginning of the file after counting.
///
/// Note: for pipeline progress estimation, [`estimate_total_lines`] is
/// preferred as it avoids a full file scan.
#[allow(dead_code)]
pub fn count_lines(path: &str) -> Result<u64> {
    let mut file = File::open(path)?;
    let count = count_lines_in_buf(&mut file)?;
    file.seek(SeekFrom::Start(0))?;
    Ok(count)
}

/// Count newlines in an already-open file using a 256KB read buffer.
fn count_lines_in_buf(file: &mut File) -> Result<u64> {
    let mut buf = [0u8; 256 * 1024];
    let mut count = 0u64;
    let mut non_empty = false;
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            if non_empty {
                count += 1; // last line without trailing newline
            }
            break;
        }
        for &b in &buf[..n] {
            if b == b'\n' {
                count += 1;
                non_empty = false;
            } else if !non_empty && b != b'\r' {
                non_empty = true;
            }
        }
    }
    Ok(count)
}

/// Print a progress bar to stderr with percentage, throughput, and ETA.
///
/// Uses `\r` to overwrite the same line. The bar width is configurable.
/// Callers should flush stderr after this to ensure immediate display.
pub fn print_progress(processed: u64, total: u64, start: &std::time::Instant, bar_width: usize) {
    let pct = if total > 0 { processed as f64 / total as f64 } else { 0.0 };
    let filled = (bar_width as f64 * pct) as usize;
    let bar = "\u{2588}".repeat(filled) + &"\u{2591}".repeat(bar_width.saturating_sub(filled));
    let elapsed = start.elapsed().as_secs_f64();
    let rate = if elapsed > 0.0 { processed as f64 / elapsed } else { 0.0 };
    let eta = if rate > 0.0 { (total.saturating_sub(processed)) as f64 / rate } else { 0.0 };

    eprint!(
        "\r    [{bar}] {pct:5.1}%  {processed}/{total}  {rate:.0} seq/s  ETA: {eta:.0}s",
        pct = pct * 100.0,
        processed = processed,
        total = total,
        rate = rate,
        eta = eta,
    );
    std::io::stderr().flush().ok();
}
