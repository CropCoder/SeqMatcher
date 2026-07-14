//! # SeqMatcher
//!
//! High-performance DNA sequence primer matching and library variant counting.
//!
//! This crate provides both a CLI tool and a library API for processing
//! large sequence files with parallel primer matching and variant counting.
//!
//! ## Library usage
//!
//! ```no_run
//! use seq_matcher::{io, matcher};
//!
//! let primer_data = io::load_primers("primers.csv")?;
//! let lib = io::load_library("library.csv", "sequence_col")?;
//! let ac_data = matcher::build_variant_ac(&lib.variants);
//! # Ok::<(), anyhow::Error>(())
//! ```
//!
//! ## CLI usage
//!
//! ```text
//! seq_matcher \
//!   --primer-csv primers.csv \
//!   --library-csv library.csv \
//!   --seq a_11:data/11_seq.txt \
//!   --output-dir output
//! ```

pub mod bio;
pub mod cli;
pub mod io;
pub mod matcher;
pub mod pipeline;
pub mod progress;
pub mod report;
pub mod types;

// Re-export commonly used types for library consumers
pub use types::{AcData, Primer, ThreadResult, Variant};
