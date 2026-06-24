mod bio;
mod io;
mod types;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;

use io::LibraryData;
use types::{Primer, ThreadResult};

#[derive(Parser, Debug)]
#[command(
    name = "seq_matcher",
    version,
    about = "高性能序列引物匹配与文库变体计数工具",
    long_about = "并行化处理海量 DNA 序列，匹配引物并统计文库变体出现频次。"
)]
struct Args {
    /// 引物 CSV 文件路径 (列: id, forward_seq, reverse_seq)
    #[arg(short = 'p', long)]
    primer_csv: String,

    /// 文库 CSV 文件路径
    #[arg(short = 'l', long)]
    library_csv: String,

    /// 文库 CSV 中序列所在列名
    #[arg(long, default_value = "single_degenerate_library_expanded_reference")]
    library_seq_col: String,

    /// 序列文件: 格式为 LABEL:PATH (如 a_11:data/11_seq.txt)
    #[arg(short = 's', long = "seq", value_parser = parse_seq_arg)]
    seq_files: Vec<SeqInput>,

    /// 输出目录
    #[arg(short = 'o', long, default_value = "output")]
    output_dir: String,

    /// 并行处理块大小 (条/批)
    #[arg(short = 'c', long, default_value = "10000")]
    chunk_size: usize,

    /// 线程数 (默认使用全部 CPU 核心)
    #[arg(short = 't', long)]
    threads: Option<usize>,
}

#[derive(Debug, Clone)]
struct SeqInput {
    label: String,
    path: String,
}

fn parse_seq_arg(s: &str) -> Result<SeqInput, String> {
    let parts: Vec<&str> = s.splitn(2, ':').collect();
    if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
        return Err(format!(
            "无效的序列参数 '{}'，期望格式 LABEL:PATH (如 a_11:data/11_seq.txt)",
            s
        ));
    }
    Ok(SeqInput {
        label: parts[0].to_string(),
        path: parts[1].to_string(),
    })
}

fn main() -> Result<()> {
    let args = Args::parse();

    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .context("无法设置线程池")?;
    }

    let output_dir = PathBuf::from(&args.output_dir);
    std::fs::create_dir_all(&output_dir)
        .with_context(|| format!("无法创建输出目录: {}", args.output_dir))?;

    eprintln!("  Loading primers from: {}", args.primer_csv);
    let primers = Arc::new(io::load_primers(&args.primer_csv)?);
    eprintln!("  Loaded {} primers", primers.len());

    eprintln!("  Loading library from: {}", args.library_csv);
    let lib = Arc::new(io::load_library(
        &args.library_csv,
        &args.library_seq_col,
    )?);
    eprintln!("  Loaded {} library variants", lib.variants.len());

    for seq_input in &args.seq_files {
        eprintln!(
            "  Processing: {} -> {}",
            seq_input.path, seq_input.label
        );
        let result = process_sequences(
            &primers,
            &lib,
            &seq_input.path,
            args.chunk_size,
        )?;

        let primer_output = output_dir.join(format!(
            "{}_seq_matched_primers_count.csv",
            seq_input.label
        ));
        let variant_output = output_dir.join(format!(
            "{}_seq_matched_library_variant_count.csv",
            seq_input.label
        ));

        io::write_primer_counts(&primer_output, &primers, &result, &seq_input.label)?;
        io::write_variant_counts(
            &variant_output,
            &lib,
            &primers,
            &result,
            &seq_input.label,
        )?;

        eprintln!(
            "  Written: {}, {}",
            primer_output.display(),
            variant_output.display()
        );
    }

    eprintln!("  All done.");
    Ok(())
}

fn process_sequences(
    primers: &[Primer],
    lib: &LibraryData,
    seq_path: &str,
    chunk_size: usize,
) -> Result<ThreadResult> {
    let file =
        File::open(seq_path).with_context(|| format!("无法打开序列文件: {}", seq_path))?;
    let reader = BufReader::new(file);

    let mut global_result = ThreadResult::default();
    let mut total = 0usize;

    let lines: Vec<String> = reader
        .lines()
        .filter_map(|l| l.ok())
        .map(|l| l.trim().to_owned())
        .filter(|l| !l.is_empty())
        .collect();

    for chunk in lines.chunks(chunk_size) {
        let chunk: Vec<&str> = chunk.iter().map(|s| s.as_str()).collect();

        let chunk_result = chunk
            .par_iter()
            .fold(ThreadResult::default, |mut local, seq| {
                let seq_upper = seq.to_uppercase();
                let matched = primers.iter().find_map(|primer| {
                    if bio::check_primer_match(&seq_upper, primer) {
                        Some(primer.id.clone())
                    } else {
                        None
                    }
                });

                if let Some(p_id) = matched {
                    *local.primer_counts.entry(p_id.clone()).or_default() += 1;
                    let var_map = local.variant_counts.entry(p_id).or_default();
                    for (idx, var) in lib.variants.iter().enumerate() {
                        if var.raw.is_empty() {
                            continue;
                        }
                        if seq_upper.contains(&var.raw) || seq_upper.contains(&var.rc) {
                            *var_map.entry(idx).or_default() += 1;
                        }
                    }
                }
                local
            })
            .reduce(ThreadResult::default, |mut a, b| {
                a.merge(&b);
                a
            });

        global_result.merge(&chunk_result);
        total += chunk.len();
        eprintln!("    processed {} sequences", total);
    }

    Ok(global_result)
}
