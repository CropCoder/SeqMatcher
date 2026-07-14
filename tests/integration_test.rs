//! Integration tests for the SeqMatcher pipeline.
//!
//! Tests the complete workflow: CSV loading → AC building → chunk processing → CSV writing.

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use seq_matcher::io;
use seq_matcher::matcher;
use seq_matcher::types::{Primer, ThreadResult, Variant};

/// Helper: create a temporary CSV file with the given content.
fn create_temp_csv(dir: &std::path::Path, name: &str, content: &str) -> std::path::PathBuf {
    let path = dir.join(name);
    let mut f = File::create(&path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
    path
}

/// Helper: create a temporary sequence file (reserved for future end-to-end tests).
#[allow(dead_code)]
fn create_temp_seq(dir: &std::path::Path, name: &str, sequences: &[&str]) -> std::path::PathBuf {
    let path = dir.join(name);
    let mut f = File::create(&path).unwrap();
    for seq in sequences {
        writeln!(f, "{}", seq).unwrap();
    }
    path
}

#[test]
fn test_load_primers_basic() {
    let dir = tempfile::tempdir().unwrap();
    let csv = create_temp_csv(
        dir.path(),
        "primers.csv",
        "id,forward,reverse\nP1,ATCG,GCTA\nP2,GGGG,CCCC\n",
    );

    let data = io::load_primers(csv.to_str().unwrap()).unwrap();
    assert_eq!(data.primers.len(), 2);
    assert_eq!(data.primers[0].id, "P1");
    assert_eq!(data.primers[0].f, "ATCG");
    assert_eq!(data.primers[0].r, "GCTA");
    assert_eq!(data.primers[1].id, "P2");
}

#[test]
fn test_load_primers_empty_id_detected() {
    // Note: empty IDs are now validated in pipeline::validate_inputs, not in load_primers
    // load_primers will accept them; validation happens later
    let dir = tempfile::tempdir().unwrap();
    let csv = create_temp_csv(dir.path(), "primers.csv", "id,forward,reverse\n,ATCG,GCTA\n");
    let data = io::load_primers(csv.to_str().unwrap()).unwrap();
    assert_eq!(data.primers[0].id, "");
}

#[test]
fn test_load_library_basic() {
    let dir = tempfile::tempdir().unwrap();
    let csv =
        create_temp_csv(dir.path(), "library.csv", "var_id,sequence\nV1,ATCGATCG\nV2,GGGGCCCC\n");

    let data = io::load_library(csv.to_str().unwrap(), "sequence").unwrap();
    assert_eq!(data.variants.len(), 2);
    assert_eq!(data.variants[0].raw, "ATCGATCG");
    assert_eq!(data.variants[1].raw, "GGGGCCCC");
}

#[test]
fn test_load_library_empty_sequence() {
    let dir = tempfile::tempdir().unwrap();
    let csv = create_temp_csv(dir.path(), "library.csv", "var_id,sequence\nV1,\nV2,ATCG\n");

    let data = io::load_library(csv.to_str().unwrap(), "sequence").unwrap();
    assert_eq!(data.variants.len(), 2);
    assert!(data.variants[0].raw.is_empty());
    assert!(!data.variants[1].raw.is_empty());
}

#[test]
fn test_build_variant_ac_empty_variants() {
    let variants = vec![Variant::new(""), Variant::new("ATCG"), Variant::new("")];
    let ac_data = matcher::build_variant_ac(&variants);

    assert_eq!(ac_data.empty_variants.len(), 2);
    assert_eq!(ac_data.empty_variants, vec![0, 2]);
    // One non-empty variant → 1 raw + 1 RC = 2 patterns
    assert_eq!(ac_data.pattern_to_variant.len(), 2);
}

#[test]
fn test_process_chunk_single_primer_match() {
    let primers = vec![Primer::new("P1".into(), "ATCG".into(), "GCTA".into())];
    let variants = vec![Variant::new("GGGG"), Variant::new("CCCC")];
    let ac_data = matcher::build_variant_ac(&variants);

    // Sequence starts with ATCG and ends with TAGC (RC of GCTA)
    let seq = "ATCGNNNNNNNNTAGC".to_string();
    let chunk = vec![seq];

    let result = matcher::process_chunk(
        &chunk,
        &primers,
        &ac_data.ac,
        &ac_data.pattern_to_variant,
        &ac_data.empty_variants,
        1,
    );

    assert_eq!(result.primer_counts.get("P1"), Some(&1));
    // No variant should be found in the sequence (GGGG and CCCC are not in seq)
    let var_map = result.variant_counts.get("P1");
    assert!(var_map.is_none() || var_map.unwrap().is_empty());
}

#[test]
fn test_process_chunk_variant_match() {
    let primers = vec![Primer::new("P1".into(), "ATCG".into(), "GCTA".into())];
    let variants = vec![Variant::new("NNNNN")];
    let ac_data = matcher::build_variant_ac(&variants);

    // Sequence contains NNNNN which is a variant
    let seq = "ATCGNNNNNTAGC".to_string();
    let chunk = vec![seq];

    let result = matcher::process_chunk(
        &chunk,
        &primers,
        &ac_data.ac,
        &ac_data.pattern_to_variant,
        &ac_data.empty_variants,
        1,
    );

    assert_eq!(result.primer_counts.get("P1"), Some(&1));
    let var_map = result.variant_counts.get("P1").unwrap();
    assert_eq!(var_map.get(&0), Some(&1)); // variant 0 matched once
}

#[test]
fn test_process_chunk_no_primer_match() {
    let primers = vec![Primer::new("P1".into(), "ATCG".into(), "GCTA".into())];
    let variants = vec![Variant::new("GGGG")];
    let ac_data = matcher::build_variant_ac(&variants);

    // Sequence does not match primer
    let seq = "GGGGNNNNCCCC".to_string();
    let chunk = vec![seq];

    let result = matcher::process_chunk(
        &chunk,
        &primers,
        &ac_data.ac,
        &ac_data.pattern_to_variant,
        &ac_data.empty_variants,
        1,
    );

    assert!(!result.primer_counts.contains_key("P1"));
}

#[test]
fn test_process_chunk_empty_variant_always_matches() {
    let primers = vec![Primer::new("P1".into(), "ATCG".into(), "GCTA".into())];
    let variants = vec![Variant::new("")]; // empty variant — always counted
    let ac_data = matcher::build_variant_ac(&variants);

    let seq = "ATCGNNNNTAGC".to_string();
    let chunk = vec![seq];

    let result = matcher::process_chunk(
        &chunk,
        &primers,
        &ac_data.ac,
        &ac_data.pattern_to_variant,
        &ac_data.empty_variants,
        1,
    );

    let var_map = result.variant_counts.get("P1").unwrap();
    assert_eq!(var_map.get(&0), Some(&1));
}

#[test]
fn test_thread_result_merge() {
    let mut a = ThreadResult::with_capacity(2);
    a.primer_counts.insert("P1".into(), 10);
    let mut inner = HashMap::new();
    inner.insert(0, 5);
    a.variant_counts.insert("P1".into(), inner);

    let mut b = ThreadResult::with_capacity(2);
    b.primer_counts.insert("P1".into(), 3);
    b.primer_counts.insert("P2".into(), 7);
    let mut inner = HashMap::new();
    inner.insert(0, 2);
    b.variant_counts.insert("P1".into(), inner);

    a.merge(&b);

    assert_eq!(a.primer_counts.get("P1"), Some(&13));
    assert_eq!(a.primer_counts.get("P2"), Some(&7));
    assert_eq!(a.variant_counts.get("P1").unwrap().get(&0), Some(&7));
}

#[test]
fn test_csv_write_roundtrip_primer_counts() {
    let dir = tempfile::tempdir().unwrap();
    let csv = create_temp_csv(
        dir.path(),
        "primers.csv",
        "id,forward,reverse\nP1,ATCG,GCTA\nP2,GGGG,CCCC\n",
    );

    let primer_data = io::load_primers(csv.to_str().unwrap()).unwrap();

    let mut result = ThreadResult::with_capacity(2);
    result.primer_counts.insert("P1".into(), 100);
    result.primer_counts.insert("P2".into(), 50);

    let output = dir.path().join("output.csv");
    io::write_primer_counts(&output, &primer_data, &result, "test").unwrap();

    // Verify the output file exists and has content
    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("count_test"));
    assert!(content.contains("100"));
    assert!(content.contains("50"));
}

#[test]
fn test_csv_write_roundtrip_variant_counts() {
    let dir = tempfile::tempdir().unwrap();
    let csv = create_temp_csv(dir.path(), "library.csv", "var_id,sequence\nV1,ATCG\nV2,GGGG\n");

    let lib = io::load_library(csv.to_str().unwrap(), "sequence").unwrap();
    let primers = vec![
        Primer::new("P1".into(), "AAAA".into(), "TTTT".into()),
        Primer::new("P2".into(), "CCCC".into(), "GGGG".into()),
    ];

    let mut result = ThreadResult::with_capacity(2);
    let mut inner1 = HashMap::new();
    inner1.insert(0, 10);
    inner1.insert(1, 20);
    result.variant_counts.insert("P1".into(), inner1);

    let output = dir.path().join("output_variants.csv");
    io::write_variant_counts(&output, &lib, &primers, &result, "test").unwrap();

    let content = std::fs::read_to_string(&output).unwrap();
    assert!(content.contains("P1_test"));
    assert!(content.contains("P2_test"));
    assert!(content.contains("10"));
    assert!(content.contains("20"));
}
