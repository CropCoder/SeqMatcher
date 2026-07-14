// SeqMatcher: High-performance DNA sequence primer matching & library variant counting.
//
// Parallelized processing of large sequence files using Rayon data parallelism
// and Aho-Corasick multi-pattern matching for O(L+M) variant detection.
//
// Usage:
//   seq_matcher -p primers.csv -l library.csv -s a_11:data/11_seq.txt

mod bio;
mod cli;
mod io;
mod matcher;
mod pipeline;
mod progress;
mod report;
mod types;

use clap::Parser;

fn main() -> anyhow::Result<()> {
    // Initialize structured logging
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .with_target(false)
        .init();

    let args = cli::Args::parse();
    pipeline::run(args)
}
