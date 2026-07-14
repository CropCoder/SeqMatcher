//! Benchmark: reverse complement computation performance.
//!
//! Compares the compile-time LUT approach against a naive match-based implementation.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use seq_matcher::bio;

/// Naive reverse complement using match statements (for comparison).
fn reverse_complement_naive(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'R' => 'Y',
            'Y' => 'R',
            'K' => 'M',
            'M' => 'K',
            'S' => 'S',
            'W' => 'W',
            'B' => 'V',
            'D' => 'H',
            'H' => 'D',
            'V' => 'B',
            'N' => 'N',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            'r' => 'y',
            'y' => 'r',
            'k' => 'm',
            'm' => 'k',
            's' => 's',
            'w' => 'w',
            'b' => 'v',
            'd' => 'h',
            'h' => 'd',
            'v' => 'b',
            'n' => 'n',
            other => other,
        })
        .collect()
}

fn bench_reverse_complement(c: &mut Criterion) {
    let mut group = c.benchmark_group("reverse_complement");
    group.sample_size(100);

    // Short sequence (typical primer length)
    let short = "ATCGRYKMSWBDHVNATCG";
    // Medium sequence (typical variant length)
    let medium = "ATCGRYKMSWBDHVN".repeat(5);
    // Long sequence (typical read length)
    let long = "ATCGRYKMSWBDHVN".repeat(50);

    group
        .bench_function("lut_short_20bp", |b| b.iter(|| bio::reverse_complement(black_box(short))));
    group.bench_function("naive_short_20bp", |b| {
        b.iter(|| reverse_complement_naive(black_box(short)))
    });

    group.bench_function("lut_medium_75bp", |b| {
        b.iter(|| bio::reverse_complement(black_box(&medium)))
    });
    group.bench_function("naive_medium_75bp", |b| {
        b.iter(|| reverse_complement_naive(black_box(&medium)))
    });

    group
        .bench_function("lut_long_750bp", |b| b.iter(|| bio::reverse_complement(black_box(&long))));
    group.bench_function("naive_long_750bp", |b| {
        b.iter(|| reverse_complement_naive(black_box(&long)))
    });

    group.finish();
}

criterion_group!(benches, bench_reverse_complement);
criterion_main!(benches);
