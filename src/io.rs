//! CSV I/O operations: loading primer and library data, writing result tables.
//!
//! Handles reading primer/libraries CSV files and writing matched count
//! output CSVs, using the `csv` crate for correct field escaping.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use csv::{Reader, StringRecord, WriterBuilder};

use crate::types::{Primer, ThreadResult, Variant};

// ---------------------------------------------------------------------------
// Loaded data containers
// ---------------------------------------------------------------------------

/// Loaded primer data: typed primer objects + original CSV records for output.
pub struct PrimerData {
    pub primers: Vec<Primer>,
    pub headers: Vec<String>,
    pub records: Vec<StringRecord>,
}

/// Loaded library data: typed variant objects + original CSV records for output.
pub struct LibraryData {
    pub variants: Vec<Variant>,
    pub headers: Vec<String>,
    pub records: Vec<StringRecord>,
}

// ---------------------------------------------------------------------------
// Loaders
// ---------------------------------------------------------------------------

/// Load primer data from a CSV file.
///
/// The CSV must have at least 3 columns: ID, forward sequence, reverse sequence.
/// Additional columns are preserved and passed through to output.
pub fn load_primers(path: &str) -> Result<PrimerData> {
    let mut reader =
        Reader::from_path(path).with_context(|| format!("Failed to read primer CSV: {}", path))?;

    let headers: Vec<String> = reader.headers()?.iter().map(|h| h.to_string()).collect();

    let mut primers = Vec::new();
    let mut records = Vec::new();
    for (i, record) in reader.records().enumerate() {
        let record = record?;
        if record.len() < 3 {
            anyhow::bail!(
                "Primer CSV row {} has insufficient fields: need at least 3 (id, forward, reverse), got {}",
                i + 1,
                record.len(),
            );
        }
        let id = record[0].trim().to_string();
        let f = record[1].trim().to_string();
        let r = record[2].trim().to_string();
        primers.push(Primer::new(id, f, r));
        records.push(record);
    }
    Ok(PrimerData { primers, headers, records })
}

/// Load library variant data from a CSV file.
///
/// The sequence column is identified by name (`seq_col_name`). Empty sequences
/// produce empty variants (always counted as a match, Python-compatible behavior).
pub fn load_library(path: &str, seq_col_name: &str) -> Result<LibraryData> {
    let mut reader =
        Reader::from_path(path).with_context(|| format!("Failed to read library CSV: {}", path))?;

    let headers: Vec<String> = reader.headers()?.iter().map(|h| h.to_string()).collect();

    let seq_col_index = headers.iter().position(|h| h == seq_col_name).with_context(|| {
        format!(
            "Column '{}' not found in library CSV. Available columns: {:?}",
            seq_col_name, headers
        )
    })?;

    let mut variants = Vec::new();
    let mut records = Vec::new();

    for result in reader.records() {
        let record = result?;
        let seq = record.get(seq_col_index).map(|s| s.trim()).unwrap_or("");
        if seq.is_empty() {
            variants.push(Variant { raw: String::new(), rc: String::new() });
        } else {
            variants.push(Variant::new(seq));
        }
        records.push(record);
    }

    Ok(LibraryData { variants, headers, records })
}

// ---------------------------------------------------------------------------
// Writers
// ---------------------------------------------------------------------------

/// Write primer match counts to a CSV file.
///
/// Output includes all original columns plus a `count_{suffix}` column
/// with the number of matched sequences per primer.
pub fn write_primer_counts(
    path: &Path,
    primer_data: &PrimerData,
    result: &ThreadResult,
    suffix: &str,
) -> Result<()> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    let mut csv_writer = WriterBuilder::new().from_writer(writer);

    // Header: original columns + count_{suffix}
    let mut header = primer_data.headers.clone();
    header.push(format!("count_{}", suffix));
    csv_writer.write_record(&header).context("Failed to write primer counts header")?;

    // Data: original row + count
    for (i, record) in primer_data.records.iter().enumerate() {
        let primer = &primer_data.primers[i];
        let count = result.primer_counts.get(&primer.id).copied().unwrap_or(0);

        let mut row: Vec<String> = record.iter().map(|f| f.to_string()).collect();
        row.push(count.to_string());
        csv_writer.write_record(&row).context("Failed to write primer counts row")?;
    }

    csv_writer.flush().context("Failed to flush primer counts CSV")?;
    Ok(())
}

/// Write variant match counts to a CSV file.
///
/// Output includes all original columns plus one column per primer
/// (`{primer_id}_{suffix}`) with the number of times each variant was
/// matched under that primer.
///
/// Uses pre-built column vectors for O(1) lookup instead of O(P*V)
/// HashMap accesses.
pub fn write_variant_counts(
    path: &Path,
    lib: &LibraryData,
    primers: &[Primer],
    result: &ThreadResult,
    suffix: &str,
) -> Result<()> {
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    let mut csv_writer = WriterBuilder::new().from_writer(writer);

    // Pre-build column data: Vec<Vec<usize>> where col_data[p][v] = count
    // This converts O(P*V) HashMap lookups to O(1) array indexing.
    let col_data: Vec<Vec<usize>> = primers
        .iter()
        .map(|p| {
            let map = result.variant_counts.get(&p.id);
            (0..lib.variants.len())
                .map(|vi| map.and_then(|m| m.get(&vi)).copied().unwrap_or(0))
                .collect()
        })
        .collect();

    // Header: original columns + per-primer count columns
    let mut header = lib.headers.clone();
    for primer in primers {
        header.push(format!("{}_{}", primer.id, suffix));
    }
    csv_writer.write_record(&header).context("Failed to write variant counts header")?;

    // Data rows
    for (idx, record) in lib.records.iter().enumerate() {
        let mut row: Vec<String> = record.iter().map(|f| f.to_string()).collect();
        for col in &col_data {
            row.push(col[idx].to_string());
        }
        csv_writer.write_record(&row).context("Failed to write variant counts row")?;
    }

    csv_writer.flush().context("Failed to flush variant counts CSV")?;
    Ok(())
}
