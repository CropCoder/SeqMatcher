//! Pipeline orchestration: ties together CLI, I/O, matching, progress,
//! and report generation into the main processing workflow.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result};
use tracing::info;

use crate::cli::Args;
use crate::matcher;
use crate::progress;
use crate::types::ThreadResult;

/// Run the full matching pipeline from parsed CLI arguments.
///
/// This is the main entry point called from `main()`. It handles:
/// 1. Thread pool configuration
/// 2. Input loading and validation
/// 3. Aho-Corasick automaton construction
/// 4. Per-file sequence processing
/// 5. Output CSV writing
/// 6. Run summary report generation
pub fn run(args: Args) -> Result<()> {
    // --- Thread pool configuration ---
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .context("Failed to configure thread pool")?;
    }

    // --- Output directory ---
    let output_dir = PathBuf::from(&args.output_dir);
    std::fs::create_dir_all(&output_dir)
        .with_context(|| format!("Failed to create output directory: {}", args.output_dir))?;

    // --- Load inputs ---
    info!("Loading primers from: {}", args.primer_csv);
    let primer_data = Arc::new(crate::io::load_primers(&args.primer_csv)?);
    info!("Loaded {} primers", primer_data.primers.len());

    info!("Loading library from: {}", args.library_csv);
    let lib = Arc::new(crate::io::load_library(&args.library_csv, &args.library_seq_col)?);
    info!("Loaded {} library variants", lib.variants.len());

    // --- Validate inputs ---
    validate_inputs(&primer_data.primers, &lib.variants, &args.seq_files)?;

    // --- Build Aho-Corasick automaton ---
    // Destructure AcData so each field can be independently wrapped in Arc
    let ac_data = matcher::build_variant_ac(&lib.variants);
    let num_patterns = ac_data.pattern_to_variant.len();
    let empty_count = ac_data.empty_variants.len();
    let ac = Arc::new(ac_data.ac);
    let pattern_to_variant = Arc::new(ac_data.pattern_to_variant);
    let empty_variants = Arc::new(ac_data.empty_variants);

    if empty_count > 0 {
        info!("Note: {} empty variant(s) — always counted per match (Python compat)", empty_count);
    }

    // --- Dry-run: stop here ---
    if args.dry_run {
        print_dry_run_summary(
            &primer_data.primers,
            &lib.variants,
            num_patterns,
            empty_count,
            &args.seq_files,
        );
        return Ok(());
    }

    // --- Process each sequence file ---
    for seq_input in &args.seq_files {
        if !args.quiet {
            info!("Processing: {} -> {}", seq_input.path, seq_input.label);
        }

        let result = process_sequences(
            &primer_data.primers,
            lib.variants.len(),
            &ac,
            &pattern_to_variant,
            &empty_variants,
            &seq_input.path,
            args.chunk_size,
            primer_data.primers.len(),
            args.quiet,
        )?;

        // --- Write outputs ---
        let primer_output =
            output_dir.join(format!("{}_seq_matched_primers_count.csv", seq_input.label));
        let variant_output =
            output_dir.join(format!("{}_seq_matched_library_variant_count.csv", seq_input.label));

        crate::io::write_primer_counts(&primer_output, &primer_data, &result, &seq_input.label)?;
        crate::io::write_variant_counts(
            &variant_output,
            &lib,
            &primer_data.primers,
            &result,
            &seq_input.label,
        )?;

        if !args.quiet {
            info!("Written: {}, {}", primer_output.display(), variant_output.display());
        }
    }

    // --- Run summary report ---
    crate::report::write_run_summary(
        &output_dir,
        &args,
        &primer_data.primers,
        lib.variants.len(),
        num_patterns,
    )?;

    info!("All done.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Input validation
// ---------------------------------------------------------------------------

/// Validate input data before processing.
///
/// Checks:
/// - No duplicate or empty primer IDs
/// - Not all variants are empty
/// - All sequence files exist and are readable
fn validate_inputs(
    primers: &[crate::types::Primer],
    variants: &[crate::types::Variant],
    seq_files: &[crate::cli::SeqInput],
) -> Result<()> {
    // Check for duplicate or empty primer IDs
    let mut seen_ids = std::collections::HashSet::new();
    for p in primers {
        if p.id.is_empty() {
            anyhow::bail!("Primer has an empty ID");
        }
        if !seen_ids.insert(&p.id) {
            anyhow::bail!("Duplicate primer ID: {}", p.id);
        }
    }

    // Check that at least some variants are non-empty
    let non_empty_count = variants.iter().filter(|v| !v.raw.is_empty()).count();
    if variants.is_empty() {
        anyhow::bail!("Library contains no variants");
    }
    if non_empty_count == 0 && !variants.is_empty() {
        info!("Warning: all library variants are empty sequences");
    }

    // Check sequence file existence
    for s in seq_files {
        let metadata = std::fs::metadata(&s.path)
            .with_context(|| format!("Sequence file not found: {}", s.path))?;
        if metadata.len() == 0 {
            anyhow::bail!("Sequence file is empty: {}", s.path);
        }
        if !metadata.is_file() {
            anyhow::bail!("Sequence path is not a regular file: {}", s.path);
        }
    }

    Ok(())
}

/// Print a validation summary for dry-run mode.
fn print_dry_run_summary(
    primers: &[crate::types::Primer],
    variants: &[crate::types::Variant],
    num_patterns: usize,
    empty_count: usize,
    seq_files: &[crate::cli::SeqInput],
) {
    info!("=== Dry-run validation summary ===");
    info!("  Primers:           {}", primers.len());
    info!("  Library variants:  {}", variants.len());
    info!("  Non-empty variants: {}", variants.iter().filter(|v| !v.raw.is_empty()).count());
    info!("  AC patterns:       {}", num_patterns);
    info!("  Empty variants:    {}", empty_count);
    info!("  Sequence files:    {}", seq_files.len());
    for s in seq_files {
        let size = std::fs::metadata(&s.path).map(|m| m.len()).unwrap_or(0);
        info!("    {} -> {} ({:.2} MB)", s.path, s.label, size as f64 / 1_048_576.0);
    }
    info!("  All inputs validated successfully.");
}

// ---------------------------------------------------------------------------
// Sequence processing
// ---------------------------------------------------------------------------

/// Process a single sequence file: stream chunks, parallel match, aggregate.
///
/// Uses file-size-based line estimation to avoid a full pre-scan.
/// Progress is reported to stderr in real-time (unless quiet mode).
#[allow(clippy::too_many_arguments)]
fn process_sequences(
    primers: &[crate::types::Primer],
    num_variants: usize,
    ac: &aho_corasick::AhoCorasick,
    pattern_to_variant: &[usize],
    empty_variants: &[usize],
    seq_path: &str,
    chunk_size: usize,
    num_primers: usize,
    quiet: bool,
) -> Result<ThreadResult> {
    // Estimate total from file size (no pre-scan)
    let (total, file_size) = progress::estimate_total_lines(seq_path, 1000)?;

    if total == 0 {
        info!("Sequence file is empty, skipping");
        return Ok(ThreadResult::with_capacity(num_primers));
    }

    if !quiet {
        info!(
            "  Total: ~{total} (est.) | {size:.1} MB | chunk: {chunk_size} | \
             primers: {p_len} | variants: {v_len} | AC patterns: {pat}",
            total = total,
            size = file_size as f64 / 1_048_576.0,
            chunk_size = chunk_size,
            p_len = primers.len(),
            v_len = num_variants,
            pat = pattern_to_variant.len(),
        );
    }

    let file = File::open(seq_path)
        .with_context(|| format!("Failed to open sequence file: {seq_path}"))?;
    let reader = BufReader::new(file);

    let mut global_result = ThreadResult::with_capacity(num_primers);
    let mut processed = 0u64;
    let bar_width = 40usize;
    let start = Instant::now();

    // Streaming chunk accumulation
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
            let chunk_result = matcher::process_chunk(
                &chunk,
                primers,
                ac,
                pattern_to_variant,
                empty_variants,
                num_primers,
            );
            global_result.merge(&chunk_result);
            processed += chunk.len() as u64;
            if !quiet {
                progress::print_progress(processed, total, &start, bar_width);
            }
            chunk.clear();
        }
    }

    // Final partial chunk
    if !chunk.is_empty() {
        let chunk_result = matcher::process_chunk(
            &chunk,
            primers,
            ac,
            pattern_to_variant,
            empty_variants,
            num_primers,
        );
        global_result.merge(&chunk_result);
        processed += chunk.len() as u64;
        if !quiet {
            progress::print_progress(processed, total, &start, bar_width);
        }
    }

    let elapsed = start.elapsed().as_secs_f64();
    info!(
        "  Complete: {} sequences in {:.1}s, {:.0} seq/s",
        processed,
        elapsed,
        processed as f64 / elapsed.max(0.001),
    );

    Ok(global_result)
}
