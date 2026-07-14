//! Benchmark: Aho-Corasick automaton construction and matching performance
//! across varying numbers of variant patterns.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use seq_matcher::matcher;
use seq_matcher::types::Variant;

/// Generate synthetic variant sequences of the given length and count.
fn gen_variants(count: usize, seq_len: usize) -> Vec<Variant> {
    let bases = ['A', 'T', 'C', 'G'];
    (0..count)
        .map(|i| {
            let mut s = String::with_capacity(seq_len);
            for j in 0..seq_len {
                let idx = (i * 31 + j * 7) % 4;
                s.push(bases[idx]);
            }
            Variant::new(&s)
        })
        .collect()
}

fn gen_sequences(count: usize, seq_len: usize) -> Vec<String> {
    let bases = ['A', 'T', 'C', 'G'];
    (0..count)
        .map(|i| {
            let mut s = String::with_capacity(seq_len);
            for j in 0..seq_len {
                let idx = (i * 17 + j * 13) % 4;
                s.push(bases[idx]);
            }
            s
        })
        .collect()
}

fn bench_ac_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("ac_build");
    group.sample_size(10);

    for &variant_count in &[100, 1000, 8192] {
        let variants = gen_variants(variant_count, 50);
        group.bench_function(format!("build_{}_variants", variant_count), |b| {
            b.iter(|| matcher::build_variant_ac(black_box(&variants)))
        });
    }
    group.finish();
}

fn bench_ac_match(c: &mut Criterion) {
    let mut group = c.benchmark_group("ac_match");
    group.sample_size(50);

    let variants = gen_variants(8192, 50);
    let ac_data = matcher::build_variant_ac(&variants);
    let sequences = gen_sequences(1000, 200);

    group.bench_function("match_1kseq_8192variants", |b| {
        b.iter(|| {
            for seq in &sequences {
                let matches: Vec<_> = ac_data
                    .ac
                    .find_iter(seq.as_bytes())
                    .map(|m| ac_data.pattern_to_variant[m.pattern().as_usize()])
                    .collect();
                black_box(matches);
            }
        })
    });
    group.finish();
}

fn bench_chunk_process(c: &mut Criterion) {
    let mut group = c.benchmark_group("chunk_process");

    let variants = gen_variants(1000, 50);
    let ac_data = matcher::build_variant_ac(&variants);
    let primers = vec![
        seq_matcher::types::Primer::new("P1".into(), "ATCG".into(), "GCTA".into()),
        seq_matcher::types::Primer::new("P2".into(), "GGGG".into(), "CCCC".into()),
    ];

    // Sequences that match the first primer
    let chunk: Vec<String> = (0..1000).map(|_| format!("ATCGAATTCCGGTAGC")).collect();

    group.bench_function("1k_seq_1k_variants", |b| {
        b.iter(|| {
            matcher::process_chunk(
                black_box(&chunk),
                black_box(&primers),
                black_box(&ac_data.ac),
                black_box(&ac_data.pattern_to_variant),
                black_box(&ac_data.empty_variants),
                2,
            )
        })
    });
    group.finish();
}

criterion_group!(benches, bench_ac_build, bench_ac_match, bench_chunk_process);
criterion_main!(benches);
