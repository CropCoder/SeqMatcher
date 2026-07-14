//! Core data types for sequence matching: primers, variants, thread-local
//! accumulators, and the Aho-Corasick build result.

use std::collections::HashMap;

use crate::bio;

// ---------------------------------------------------------------------------
// Primer & Variant
// ---------------------------------------------------------------------------

/// A primer pair with pre-computed reverse complements.
///
/// The forward (`f`) and reverse (`r`) sequences are stored in uppercase.
/// Reverse complements (`rc_f`, `rc_r`) are computed at construction time
/// to avoid repeated computation during matching.
#[derive(Debug, Clone)]
pub struct Primer {
    /// Unique primer identifier (from CSV first column)
    pub id: String,
    /// Forward primer sequence (uppercase)
    pub f: String,
    /// Reverse primer sequence (uppercase)
    pub r: String,
    /// Reverse-complement of the forward sequence
    pub rc_f: String,
    /// Reverse-complement of the reverse sequence
    pub rc_r: String,
}

impl Primer {
    /// Create a new primer, automatically computing reverse complements
    /// and converting sequences to uppercase.
    pub fn new(id: String, f: String, r: String) -> Self {
        let rc_f = bio::reverse_complement(&f);
        let rc_r = bio::reverse_complement(&r);
        Self { id, f: f.to_uppercase(), r: r.to_uppercase(), rc_f, rc_r }
    }
}

/// A library variant with pre-computed reverse complement.
///
/// Both `raw` and `rc` are stored in uppercase. An empty `raw` field
/// represents a "null" variant that always counts as matched (for
/// compatibility with Python reference implementation).
#[derive(Debug, Clone)]
pub struct Variant {
    /// Original variant sequence (uppercase)
    pub raw: String,
    /// Reverse-complement of the variant sequence (uppercase)
    pub rc: String,
}

impl Variant {
    /// Create a new variant from a sequence string.
    /// Automatically converts to uppercase and computes the reverse complement.
    pub fn new(seq: &str) -> Self {
        let raw = seq.to_uppercase();
        let rc = bio::reverse_complement(&raw);
        Self { raw, rc }
    }
}

// ---------------------------------------------------------------------------
// Thread-local accumulator
// ---------------------------------------------------------------------------

/// Per-thread accumulator for parallel sequence processing.
///
/// Each Rayon thread accumulates counts locally (lock-free), then results
/// are merged via [`ThreadResult::merge`] into a global accumulator.
///
/// `hit_buf` is a reusable scratch buffer used transiently during fold
/// for AC match deduplication — it is **never** merged, only cleared
/// and reused per-sequence.
#[derive(Clone, Default)]
pub struct ThreadResult {
    /// `primer_id → match_count`
    pub primer_counts: HashMap<String, usize>,
    /// `primer_id → (variant_index → match_count)`
    pub variant_counts: HashMap<String, HashMap<usize, usize>>,
    /// Reusable buffer for AC match dedup (thread-local, never merged)
    pub hit_buf: Vec<usize>,
}

impl ThreadResult {
    /// Create a new accumulator with pre-allocated capacity for the given
    /// number of primers.
    pub fn with_capacity(num_primers: usize) -> Self {
        Self {
            primer_counts: HashMap::with_capacity(num_primers),
            variant_counts: HashMap::with_capacity(num_primers),
            hit_buf: Vec::new(),
        }
    }

    /// Merge another thread's results into this one.
    ///
    /// Adds counts from `other` into `self`. The `hit_buf` is not merged
    /// (it is a thread-local scratch buffer with no merge semantics).
    pub fn merge(&mut self, other: &ThreadResult) {
        for (k, v) in &other.primer_counts {
            *self.primer_counts.entry(k.clone()).or_default() += v;
        }
        for (p_id, vars) in &other.variant_counts {
            let entry = self.variant_counts.entry(p_id.clone()).or_default();
            for (v_idx, count) in vars {
                *entry.entry(*v_idx).or_default() += count;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Aho-Corasick build result
// ---------------------------------------------------------------------------

/// Result of building the Aho-Corasick automaton from library variants.
pub struct AcData {
    /// The compiled Aho-Corasick automaton
    pub ac: aho_corasick::AhoCorasick,
    /// Maps AC pattern index → variant index
    pub pattern_to_variant: Vec<usize>,
    /// Indices of variants with empty sequences (always counted as a match)
    pub empty_variants: Vec<usize>,
}
