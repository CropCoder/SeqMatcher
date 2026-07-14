//! Run summary report generation.
//!
//! Produces machine-readable (JSON) and human-readable (TXT) reports
//! documenting the complete run configuration, input parameters, system
//! information, and aggregate statistics — essential for academic
//! reproducibility and pipeline auditing.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::cli::Args;
use crate::types::Primer;

/// Write both `run_summary.json` and `run_summary.txt` to the output directory.
pub fn write_run_summary(
    output_dir: &Path,
    args: &Args,
    primers: &[Primer],
    num_variants: usize,
    num_patterns: usize,
) -> Result<()> {
    let timestamp = timestamp_iso8601();

    // JSON (machine-readable)
    let json_path = output_dir.join("run_summary.json");
    let json = build_json_summary(args, primers, num_variants, num_patterns, &timestamp);
    let mut f = BufWriter::new(File::create(&json_path)?);
    writeln!(f, "{}", json)?;
    f.flush()?;

    // TXT (human-readable)
    let txt_path = output_dir.join("run_summary.txt");
    let txt = build_text_summary(args, primers, num_variants, num_patterns, &timestamp);
    let mut f = BufWriter::new(File::create(&txt_path)?);
    write!(f, "{}", txt)?;
    f.flush()?;

    Ok(())
}

// ---------------------------------------------------------------------------
// Timestamp
// ---------------------------------------------------------------------------

/// Generate an ISO 8601 timestamp string (basic format, UTC).
fn timestamp_iso8601() -> String {
    let duration =
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap_or_default();
    let secs = duration.as_secs();
    // Convert UNIX timestamp to a readable UTC datetime string
    let days_since_epoch = secs / 86400;
    let time_of_day = secs % 86400;
    let hours = time_of_day / 3600;
    let minutes = (time_of_day % 3600) / 60;
    let seconds = time_of_day % 60;

    // Simple civil date calculation from days since UNIX epoch
    // This is correct for dates after 1970-01-01
    let (year, month, day) = civil_from_days(days_since_epoch as i64);

    format!("{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z", year, month, day, hours, minutes, seconds)
}

/// Convert days since 1970-01-01 to (year, month, day).
/// Uses the algorithm from Howard Hinnant's date library (public domain).
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u32; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; // [0, 399]
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = doy - (153 * mp + 2) / 5 + 1; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 }; // [1, 12]
    let y = y + if m <= 2 { 1 } else { 0 };
    (y, m, d)
}

// ---------------------------------------------------------------------------
// JSON summary
// ---------------------------------------------------------------------------

fn build_json_summary(
    args: &Args,
    primers: &[Primer],
    num_variants: usize,
    num_patterns: usize,
    timestamp: &str,
) -> String {
    let git_commit = option_env!("SEQMATCHER_GIT_HASH").unwrap_or("unknown");
    let version = env!("CARGO_PKG_VERSION");
    let cpu_cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);

    let seq_files_json: Vec<String> = args
        .seq_files
        .iter()
        .map(|s| {
            let size = std::fs::metadata(&s.path).map(|m| m.len()).unwrap_or(0);
            format!(
                r#"      {{ "label": "{}", "path": "{}", "size_bytes": {} }}"#,
                s.label, s.path, size
            )
        })
        .collect();

    format!(
        r#"{{
  "tool": "SeqMatcher",
  "version": "{version}",
  "git_commit": "{git_commit}",
  "timestamp": "{timestamp}",
  "input": {{
    "primer_csv": "{primer_csv}",
    "primer_count": {primer_count},
    "library_csv": "{library_csv}",
    "variant_count": {variant_count},
    "ac_pattern_count": {num_patterns},
    "seq_files": [
{seq_files}
    ]
  }},
  "parameters": {{
    "chunk_size": {chunk_size},
    "threads": {threads},
    "dry_run": {dry_run}
  }},
  "system": {{
    "platform": "{platform}",
    "cpu_cores": {cpu_cores}
  }}
}}"#,
        primer_csv = args.primer_csv,
        primer_count = primers.len(),
        library_csv = args.library_csv,
        variant_count = num_variants,
        seq_files = seq_files_json.join(",\n"),
        chunk_size = args.chunk_size,
        threads = args.threads.map_or("auto".to_string(), |n| n.to_string()),
        dry_run = args.dry_run,
        platform = std::env::consts::OS,
    )
}

// ---------------------------------------------------------------------------
// Human-readable text summary
// ---------------------------------------------------------------------------

fn build_text_summary(
    args: &Args,
    primers: &[Primer],
    num_variants: usize,
    num_patterns: usize,
    timestamp: &str,
) -> String {
    let cpu_cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);

    let mut s = String::new();

    s.push_str("SeqMatcher Run Summary\n");
    s.push_str("======================\n\n");
    s.push_str(&format!("  Tool version: {}\n", env!("CARGO_PKG_VERSION")));
    s.push_str(&format!("  Timestamp:    {}\n", timestamp));
    s.push_str(&format!(
        "  Git commit:   {}\n\n",
        option_env!("SEQMATCHER_GIT_HASH").unwrap_or("unknown")
    ));

    s.push_str("Input\n");
    s.push_str("-----\n");
    s.push_str(&format!("  Primer CSV:    {} ({} primers)\n", args.primer_csv, primers.len()));
    s.push_str(&format!(
        "  Library CSV:   {} ({} variants, {} AC patterns)\n",
        args.library_csv, num_variants, num_patterns
    ));
    s.push_str("  Sequence files:\n");
    for seq in &args.seq_files {
        let size = std::fs::metadata(&seq.path).map(|m| m.len()).unwrap_or(0);
        s.push_str(&format!(
            "    {} -> {} ({:.2} MB)\n",
            seq.path,
            seq.label,
            size as f64 / 1_048_576.0
        ));
    }
    s.push('\n');

    s.push_str("Parameters\n");
    s.push_str("----------\n");
    s.push_str(&format!("  Chunk size:  {}\n", args.chunk_size));
    s.push_str(&format!(
        "  Threads:     {}\n",
        args.threads.map_or("auto".to_string(), |n| n.to_string())
    ));
    s.push_str(&format!("  Output dir:  {}\n", args.output_dir));
    s.push('\n');

    s.push_str("System\n");
    s.push_str("------\n");
    s.push_str(&format!("  Platform:  {}\n", std::env::consts::OS));
    s.push_str(&format!("  CPU cores: {}\n", cpu_cores));
    s.push('\n');

    s.push_str("Reproducibility\n");
    s.push_str("---------------\n");
    s.push_str("  To reproduce this run:\n");
    let seq_args = args
        .seq_files
        .iter()
        .map(|s| format!("{}:{}", s.label, s.path))
        .collect::<Vec<_>>()
        .join(" -s ");
    s.push_str(&format!(
        "    seq_matcher -p {} -l {} -s {} -o {} -c {}",
        args.primer_csv, args.library_csv, seq_args, args.output_dir, args.chunk_size,
    ));
    if let Some(t) = args.threads {
        s.push_str(&format!(" -t {}", t));
    }
    s.push('\n');

    s
}
