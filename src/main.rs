mod bio;
mod io;
mod types;

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::Arc;

use aho_corasick::AhoCorasick;
use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;

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
    #[arg(short = 'c', long, default_value = "100000")]
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
    let primer_data = Arc::new(io::load_primers(&args.primer_csv)?);
    eprintln!("  Loaded {} primers", primer_data.primers.len());

    eprintln!("  Loading library from: {}", args.library_csv);
    let lib = Arc::new(io::load_library(&args.library_csv, &args.library_seq_col)?);
    eprintln!("  Loaded {} library variants", lib.variants.len());

    let (ac, pattern_to_variant, empty_variants) = build_variant_ac(&lib.variants);
    let ac = Arc::new(ac);
    if !empty_variants.is_empty() {
        eprintln!(
            "  Note: {} empty variant(s) — always counted per match (Python compat)",
            empty_variants.len()
        );
    }

    let num_primers = primer_data.primers.len();

    for seq_input in &args.seq_files {
        eprintln!(
            "  Processing: {} -> {}",
            seq_input.path, seq_input.label
        );
        let result = process_sequences(
            &primer_data.primers,
            &lib.variants.len(),
            &ac,
            &pattern_to_variant,
            &empty_variants,
            &seq_input.path,
            args.chunk_size,
            num_primers,
        )?;

        let primer_output = output_dir.join(format!(
            "{}_seq_matched_primers_count.csv",
            seq_input.label
        ));
        let variant_output = output_dir.join(format!(
            "{}_seq_matched_library_variant_count.csv",
            seq_input.label
        ));

        io::write_primer_counts(&primer_output, &primer_data, &result, &seq_input.label)?;
        io::write_variant_counts(
            &variant_output,
            &lib,
            &primer_data.primers,
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

/// Build Aho-Corasick automaton from variant patterns.
/// Returns (ac, pattern→variant mapping, empty variant indices).
fn build_variant_ac(variants: &[types::Variant]) -> (AhoCorasick, Vec<usize>, Vec<usize>) {
    let mut patterns = Vec::new();
    let mut pattern_to_variant = Vec::new();
    let mut empty_variants = Vec::new();

    for (i, var) in variants.iter().enumerate() {
        if var.raw.is_empty() {
            empty_variants.push(i);
            continue;
        }
        patterns.push(var.raw.clone());
        pattern_to_variant.push(i);
        if !var.rc.is_empty() && var.rc != var.raw {
            patterns.push(var.rc.clone());
            pattern_to_variant.push(i);
        }
    }

    let ac = if patterns.is_empty() {
        AhoCorasick::new(["__SEQMATCHER_NOOP__"]).unwrap()
    } else {
        AhoCorasick::new(&patterns).expect("failed to build Aho-Corasick automaton")
    };

    (ac, pattern_to_variant, empty_variants)
}

/// Count newlines in a file for progress bar total.
fn count_lines(path: &str) -> Result<u64> {
    let mut file = File::open(path)?;
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
    file.seek(SeekFrom::Start(0))?;
    Ok(count)
}

fn process_sequences(
    primers: &[Primer],
    num_variants: &usize,
    ac: &AhoCorasick,
    pattern_to_variant: &[usize],
    empty_variants: &[usize],
    seq_path: &str,
    chunk_size: usize,
    num_primers: usize,
) -> Result<ThreadResult> {
    eprint!("    正在统计总条数... ");
    std::io::stderr().flush().ok();
    let total = count_lines(seq_path)?;
    eprintln!("{total}");

    if total == 0 {
        eprintln!("    序列文件为空，跳过");
        return Ok(ThreadResult::with_capacity(num_primers));
    }

    eprintln!(
        "    总条数: {total}  |  chunk: {chunk_size}  |  引物: {p_len}  |  变体: {v_len}  |  AC patterns: {pat}",
        total = total,
        chunk_size = chunk_size,
        p_len = primers.len(),
        v_len = num_variants,
        pat = pattern_to_variant.len(),
    );

    let file = File::open(seq_path)
        .with_context(|| format!("无法打开序列文件: {seq_path}"))?;
    let reader = BufReader::new(file);

    let mut global_result = ThreadResult::with_capacity(num_primers);
    let mut processed = 0u64;
    let bar_width = 40usize;
    let start = std::time::Instant::now();

    // Streaming: read chunk_size lines, process, free, repeat
    let mut chunk: Vec<String> = Vec::with_capacity(chunk_size);

    for line_res in reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let mut s = trimmed.to_owned();
        s.make_ascii_uppercase(); // in-place, no extra allocation
        chunk.push(s);

        if chunk.len() >= chunk_size {
            let chunk_result = process_chunk(
                &chunk, primers, ac, pattern_to_variant, empty_variants, num_primers,
            );
            global_result.merge(&chunk_result);
            processed += chunk.len() as u64;
            print_progress(processed, total, &start, bar_width);
            chunk.clear(); // free the strings
        }
    }

    // Process final partial chunk
    if !chunk.is_empty() {
        let chunk_result = process_chunk(
            &chunk, primers, ac, pattern_to_variant, empty_variants, num_primers,
        );
        global_result.merge(&chunk_result);
        processed += chunk.len() as u64;
        print_progress(processed, total, &start, bar_width);
    }

    let elapsed = start.elapsed().as_secs_f64();
    eprintln!(
        "\n    完成: {total} 条序列, 耗时 {elapsed:.1}s, 速度 {rate:.0} seq/s",
        total = total,
        elapsed = elapsed,
        rate = total as f64 / elapsed.max(0.001),
    );

    Ok(global_result)
}

/// Process a single chunk with rayon parallel fold.
fn process_chunk(
    chunk: &[String],
    primers: &[Primer],
    ac: &AhoCorasick,
    pattern_to_variant: &[usize],
    empty_variants: &[usize],
    num_primers: usize,
) -> ThreadResult {
    chunk
        .par_iter()
        .fold(
            || ThreadResult::with_capacity(num_primers),
            |mut local, seq| {
                // primer search (first match wins)
                let matched = primers.iter().find_map(|primer| {
                    if bio::check_primer_match(seq, primer) {
                        Some(primer.id.clone())
                    } else {
                        None
                    }
                });

                if let Some(p_id) = matched {
                    *local.primer_counts.entry(p_id.clone()).or_default() += 1;
                    let var_map = local.variant_counts.entry(p_id).or_default();

                    // Reuse hit_buf to avoid per-sequence allocation
                    local.hit_buf.clear();
                    local.hit_buf.extend(
                        ac.find_iter(seq.as_bytes())
                            .map(|m| pattern_to_variant[m.pattern().as_usize()]),
                    );
                    local.hit_buf.sort_unstable();
                    local.hit_buf.dedup();
                    for &vi in &local.hit_buf {
                        *var_map.entry(vi).or_default() += 1;
                    }

                    // Empty variants always match
                    for &vi in empty_variants {
                        *var_map.entry(vi).or_default() += 1;
                    }
                }
                local
            },
        )
        .reduce(
            || ThreadResult::with_capacity(num_primers),
            |mut a, b| {
                a.merge(&b);
                a
            },
        )
}

fn print_progress(processed: u64, total: u64, start: &std::time::Instant, bar_width: usize) {
    let pct = processed as f64 / total as f64;
    let filled = (bar_width as f64 * pct) as usize;
    let bar = "█".repeat(filled) + &"░".repeat(bar_width - filled);
    let elapsed = start.elapsed().as_secs_f64();
    let rate = if elapsed > 0.0 { processed as f64 / elapsed } else { 0.0 };
    let eta = if rate > 0.0 { (total - processed) as f64 / rate } else { 0.0 };

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

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_variant_search(seq: &str, variants: &[types::Variant]) -> Vec<usize> {
        variants
            .iter()
            .enumerate()
            .filter(|(_, v)| {
                if v.raw.is_empty() {
                    return true;
                }
                seq.contains(&v.raw) || (!v.rc.is_empty() && seq.contains(&v.rc))
            })
            .map(|(i, _)| i)
            .collect()
    }

    #[test]
    fn test_ac_matches_naive() {
        let variants: Vec<types::Variant> = vec![
            types::Variant::new("ATCGAAAA"),
            types::Variant::new("GGGCCCCC"),
            types::Variant::new("TTTTAAAA"),
            types::Variant::new(""),
        ];

        let (ac, p2v, empty) = build_variant_ac(&variants);
        assert_eq!(empty, vec![3]);

        let seq = "NNNATCGAAAANNNGGGCCCCCNNN";

        let mut ac_hits: Vec<usize> =
            ac.find_iter(seq).map(|m| p2v[m.pattern().as_usize()]).collect();
        ac_hits.sort_unstable();
        ac_hits.dedup();
        for &vi in &empty {
            ac_hits.push(vi);
        }
        ac_hits.sort_unstable();
        ac_hits.dedup();

        let mut naive = naive_variant_search(seq, &variants);
        naive.sort_unstable();

        assert_eq!(ac_hits, naive, "Aho-Corasick must match naive contains()");
    }

    #[test]
    fn test_ac_reverse_complement_match() {
        let variants = vec![types::Variant::new("ATCG")];
        let (ac, p2v, empty) = build_variant_ac(&variants);
        assert!(empty.is_empty());

        let seq = "NNNNCGATNNNN";
        let hits: Vec<_> = ac.find_iter(seq).map(|m| p2v[m.pattern().as_usize()]).collect();
        assert!(!hits.is_empty(), "AC should find RC match");
    }

    #[test]
    fn test_count_lines() {
        use std::io::Write;
        let dir = std::env::temp_dir();
        let path = dir.join("seqmatcher_test_count.txt");
        let mut f = File::create(&path).unwrap();
        writeln!(f, "line1\nline2\nline3").unwrap();
        drop(f);
        assert_eq!(count_lines(path.to_str().unwrap()).unwrap(), 3);
        std::fs::remove_file(&path).ok();
    }
}
