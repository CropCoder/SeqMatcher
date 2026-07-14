//! Core matching engine: Aho-Corasick automaton construction and
//! parallel chunk processing with primer search and variant counting.

use aho_corasick::AhoCorasick;
use rayon::prelude::*;

use crate::bio;
use crate::types::{AcData, Primer, ThreadResult, Variant};

/// Build an Aho-Corasick automaton from all library variant sequences.
///
/// Each non-empty variant contributes its raw sequence as a pattern, and
/// its reverse complement (if different) as an additional pattern. Both
/// patterns map to the same variant index via `pattern_to_variant`.
///
/// Empty variants are collected separately — they are always counted as
/// a match for compatibility with the Python reference implementation.
///
/// # Panics
///
/// Panics if the Aho-Corasick automaton cannot be built from the patterns
/// (should not occur for valid DNA sequences).
pub fn build_variant_ac(variants: &[Variant]) -> AcData {
    let mut patterns = Vec::new();
    let mut pattern_to_variant = Vec::new();
    let mut empty_variants = Vec::new();

    for (i, var) in variants.iter().enumerate() {
        if var.raw.is_empty() {
            empty_variants.push(i);
            continue;
        }
        // Raw sequence pattern
        patterns.push(var.raw.clone());
        pattern_to_variant.push(i);
        // Reverse complement pattern (if distinct)
        if !var.rc.is_empty() && var.rc != var.raw {
            patterns.push(var.rc.clone());
            pattern_to_variant.push(i);
        }
    }

    let ac = if patterns.is_empty() {
        // AhoCorasick::new requires at least one pattern; use a sentinel
        AhoCorasick::new(["__SEQMATCHER_NOOP__"]).expect("failed to build Aho-Corasick automaton")
    } else {
        AhoCorasick::new(&patterns).expect("failed to build Aho-Corasick automaton")
    };

    AcData { ac, pattern_to_variant, empty_variants }
}

/// Process a single chunk of sequences with Rayon parallel fold.
///
/// Each sequence is tested against all primers (first-match-wins).
/// When a primer matches, the Aho-Corasick automaton scans the sequence
/// once to find all matching variants, plus any empty variants.
///
/// # Performance
///
/// - Primer search: O(num_primers) per sequence (linear scan).
///   For typical primer sets (<100), this is negligible.
///   TODO: For >500 primers, consider a prefix-trie index.
/// - Variant search: O(L + Z) via Aho-Corasick (L=seq length, Z=matches).
/// - Deduplication: sort + dedup on matched variant indices.
pub fn process_chunk(
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
                // Primer search: first match wins
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

                    // AC scan + dedup via reusable scratch buffer
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

                    // Empty variants are always counted
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Naive O(V * L) variant search for use as a correctness reference.
    fn naive_variant_search(seq: &str, variants: &[Variant]) -> Vec<usize> {
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
        let variants: Vec<Variant> = vec![
            Variant::new("ATCGAAAA"),
            Variant::new("GGGCCCCC"),
            Variant::new("TTTTAAAA"),
            Variant::new(""), // empty variant
        ];

        let ac_data = build_variant_ac(&variants);
        assert_eq!(ac_data.empty_variants, vec![3]);

        let seq = "NNNATCGAAAANNNGGGCCCCCNNN";

        let mut ac_hits: Vec<usize> = ac_data
            .ac
            .find_iter(seq)
            .map(|m| ac_data.pattern_to_variant[m.pattern().as_usize()])
            .collect();
        ac_hits.sort_unstable();
        ac_hits.dedup();
        for &vi in &ac_data.empty_variants {
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
        let variants = vec![Variant::new("ATCG")];
        let ac_data = build_variant_ac(&variants);
        assert!(ac_data.empty_variants.is_empty());

        let seq = "NNNNCGATNNNN"; // contains CGAT = RC of ATCG
        let hits: Vec<_> = ac_data
            .ac
            .find_iter(seq)
            .map(|m| ac_data.pattern_to_variant[m.pattern().as_usize()])
            .collect();
        assert!(!hits.is_empty(), "AC should find reverse-complement match");
    }

    #[test]
    fn test_build_ac_with_all_empty() {
        let variants = vec![Variant::new("")];
        let ac_data = build_variant_ac(&variants);
        assert_eq!(ac_data.empty_variants, vec![0]);
        // Should build successfully with sentinel pattern
    }
}
