use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use csv::{Reader, StringRecord};

use crate::types::{Primer, Variant, ThreadResult};

pub struct PrimerData {
    pub primers: Vec<Primer>,
    pub headers: Vec<String>,
    pub records: Vec<StringRecord>,
}

pub struct LibraryData {
    pub variants: Vec<Variant>,
    pub headers: Vec<String>,
    pub records: Vec<StringRecord>,
}

pub fn load_primers(path: &str) -> Result<PrimerData> {
    let mut reader = Reader::from_path(path)
        .with_context(|| format!("无法读取引物 CSV: {}", path))?;

    let headers: Vec<String> = reader
        .headers()?
        .iter()
        .map(|h| h.to_string())
        .collect();

    let mut primers = Vec::new();
    let mut records = Vec::new();
    for (i, record) in reader.records().enumerate() {
        let record = record?;
        if record.len() < 3 {
            anyhow::bail!(
                "引物 CSV 第 {} 行字段不足: 需要至少 3 列 (id, forward, reverse)",
                i + 1
            );
        }
        let id = record[0].trim().to_string();
        let f = record[1].trim().to_string();
        let r = record[2].trim().to_string();
        primers.push(Primer::new(id, f, r));
        records.push(record);
    }
    Ok(PrimerData {
        primers,
        headers,
        records,
    })
}

pub fn load_library(path: &str, seq_col_name: &str) -> Result<LibraryData> {
    let mut reader = Reader::from_path(path)
        .with_context(|| format!("无法读取库 CSV: {}", path))?;

    let headers: Vec<String> = reader
        .headers()?
        .iter()
        .map(|h| h.to_string())
        .collect();

    let seq_col_index = headers
        .iter()
        .position(|h| h == seq_col_name)
        .with_context(|| {
            format!(
                "库 CSV 中找不到列 '{}'，可用列: {:?}",
                seq_col_name, headers
            )
        })?;

    let mut variants = Vec::new();
    let mut records = Vec::new();

    for result in reader.records() {
        let record = result?;
        let seq = record
            .get(seq_col_index)
            .map(|s| s.trim())
            .unwrap_or("");
        if seq.is_empty() {
            variants.push(Variant {
                raw: String::new(),
                rc: String::new(),
            });
        } else {
            variants.push(Variant::new(seq));
        }
        records.push(record);
    }

    Ok(LibraryData {
        variants,
        headers,
        records,
    })
}

pub fn write_primer_counts(
    path: &Path,
    primer_data: &PrimerData,
    result: &ThreadResult,
    suffix: &str,
) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // header: all original columns + count_{suffix}
    let header_line = primer_data.headers.join(",");
    writeln!(writer, "{},count_{}", header_line, suffix)?;

    // data: original row + count
    for (i, record) in primer_data.records.iter().enumerate() {
        let primer = &primer_data.primers[i];
        let count = result.primer_counts.get(&primer.id).copied().unwrap_or(0);
        writeln!(writer, "{},{}", csv_to_line(record), count)?;
    }
    writer.flush()?;
    Ok(())
}

pub fn write_variant_counts(
    path: &Path,
    lib: &LibraryData,
    primers: &[Primer],
    result: &ThreadResult,
    suffix: &str,
) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Build headers: original + per-primer count columns
    let count_headers: Vec<String> = primers
        .iter()
        .map(|p| format!("{}_{}", p.id, suffix))
        .collect();

    // Write header row
    let mut header_line = lib.headers.join(",");
    for ch in &count_headers {
        header_line.push(',');
        header_line.push_str(ch);
    }
    writeln!(writer, "{}", header_line)?;

    // Write data rows
    for (idx, record) in lib.records.iter().enumerate() {
        write!(writer, "{}", csv_to_line(record))?;
        for primer in primers {
            let count = result
                .variant_counts
                .get(&primer.id)
                .and_then(|m| m.get(&idx))
                .copied()
                .unwrap_or(0);
            write!(writer, ",{}", count)?;
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

fn csv_to_line(record: &StringRecord) -> String {
    record.iter().fold(String::new(), |mut acc, field| {
        if !acc.is_empty() {
            acc.push(',');
        }
        if field.contains(',') || field.contains('"') || field.contains('\n') {
            acc.push('"');
            acc.push_str(&field.replace('"', "\"\""));
            acc.push('"');
        } else {
            acc.push_str(field);
        }
        acc
    })
}
