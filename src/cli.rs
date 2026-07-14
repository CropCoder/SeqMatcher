//! Command-line interface definition using clap derive.
//!
//! Defines the argument structure and custom parsers for the
//! `LABEL:PATH` sequence file input format.

use clap::Parser;

/// High-performance DNA sequence primer matching and library variant counting tool.
///
/// Processes large sequence files in parallel, matching reads against primer pairs
/// and counting library variant occurrences using Aho-Corasick multi-pattern matching.
#[derive(Parser, Debug)]
#[command(
    name = "seq_matcher",
    version,
    about = "High-performance sequence primer matching & variant counting",
    long_about = "Parallelized DNA sequence processing: match reads to primer pairs \
                  and quantify library variant frequencies using Aho-Corasick \
                  multi-pattern matching with Rayon parallelism."
)]
pub struct Args {
    /// Primer CSV file path (columns: id, forward_seq, reverse_seq, ...)
    #[arg(short = 'p', long)]
    pub primer_csv: String,

    /// Library variant CSV file path
    #[arg(short = 'l', long)]
    pub library_csv: String,

    /// Column name in the library CSV that contains the variant sequence
    #[arg(long, default_value = "single_degenerate_library_expanded_reference")]
    pub library_seq_col: String,

    /// Sequence files in LABEL:PATH format (e.g. a_11:data/11_seq.txt).
    /// Can be specified multiple times for batch processing.
    #[arg(short = 's', long = "seq", value_parser = parse_seq_arg)]
    pub seq_files: Vec<SeqInput>,

    /// Output directory for result CSV files
    #[arg(short = 'o', long, default_value = "output")]
    pub output_dir: String,

    /// Number of sequences per parallel processing chunk
    #[arg(short = 'c', long, default_value = "100000")]
    pub chunk_size: usize,

    /// Number of worker threads (default: all available CPU cores)
    #[arg(short = 't', long)]
    pub threads: Option<usize>,

    /// Validate inputs and exit without processing sequences
    #[arg(long)]
    pub dry_run: bool,

    /// Suppress non-error output
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Enable verbose debug output
    #[arg(short = 'v', long)]
    pub verbose: bool,

    /// Append timestamp to output filenames for reproducibility
    #[arg(long)]
    pub timestamp_output: bool,
}

/// A labeled sequence file input.
///
/// The label is used for output file naming and report identification.
#[derive(Debug, Clone)]
pub struct SeqInput {
    /// Short label identifying this input (e.g. "a_11")
    pub label: String,
    /// Filesystem path to the sequence file
    pub path: String,
}

/// Parse a `LABEL:PATH` formatted sequence argument.
///
/// The first colon separates the label from the path.
/// Both parts must be non-empty.
fn parse_seq_arg(s: &str) -> Result<SeqInput, String> {
    let parts: Vec<&str> = s.splitn(2, ':').collect();
    if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
        return Err(format!(
            "Invalid sequence argument '{}': expected LABEL:PATH format (e.g. a_11:data/11_seq.txt)",
            s
        ));
    }
    Ok(SeqInput { label: parts[0].to_string(), path: parts[1].to_string() })
}
