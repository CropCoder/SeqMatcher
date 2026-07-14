# SeqMatcher

High-performance multi-threaded DNA sequence primer matching and library variant counting tool.

Matches massive sequencing reads against known primer pairs, quantifies per-primer sequence coverage, and counts library variant occurrences within matched sequences. Full Rust implementation with Rayon parallel processing and Aho-Corasick multi-pattern matching.

## Installation

```bash
git clone https://github.com/CropCoder/SeqMatcher.git
cd SeqMatcher
cargo build --release
```

The binary will be at `target/release/seq_matcher`.

### Cross-compile from macOS to Linux

```bash
# Prerequisites: brew install musl-cross && rustup target add x86_64-unknown-linux-musl
./build-linux.sh
```

Produces `target/x86_64-unknown-linux-musl/release/seq_matcher` — a statically-linked ELF 64-bit binary.

## Quick Start

```bash
./target/release/seq_matcher \
  --primer-csv primers.csv \
  --library-csv library.csv \
  --seq a_11:data/11_seq.txt \
  --output-dir output
```

### Example Run

```
  Loading primers from: primers_list_all.csv
  Loaded 32 primers
  Loading library from: 80_full_library_2_12.csv
  Loaded 8192 library variants
  Total: ~1500000 (est.) | 2.7 MB | chunk: 100000 | primers: 32 | variants: 8192 | AC patterns: 16384
  [████████████████████░░░░░░░░░░░░░░░░]  55.0%  825000/1500000  45230 seq/s  ETA: 15s
  Complete: 1500000 sequences in 33.2s, 45181 seq/s
  Written: output/a_11_seq_matched_primers_count.csv, output/a_11_seq_matched_library_variant_count.csv
  All done.
```

## Input File Formats

**Primer CSV** (`--primer-csv`) — any number of columns; all original columns are preserved in output. The first three columns are treated as primer ID, forward sequence, and reverse sequence.

| primer_id | forward_seq | reverse_seq | extra_columns... |
|-----------|-------------|-------------|-------------------|
| P001      | ATCGGTACC   | GCTATAGCA   | (preserved)        |
| P002      | TGCACTGAC   | CGTACGATG   | (preserved)        |

**Library CSV** (`--library-csv`) — contains a variant sequence column whose name is configured via `--library-seq-col`.

| variant_id | single_degenerate_library_expanded_reference | extra_columns... |
|------------|---------------------------------------------|-------------------|
| V001       | ATCGNNNTCGA                                 | (preserved)        |

**Sequence files** (`--seq`) — plain text, one DNA sequence per line.

## Output Files

Two CSV files are generated per `--seq` input:

- `{LABEL}_seq_matched_primers_count.csv` — original primer table with an added `count_{LABEL}` column (number of sequences matched to each primer)
- `{LABEL}_seq_matched_library_variant_count.csv` — original library table with per-primer variant count columns (`{primer_id}_{LABEL}`)

Additionally, each run produces:

- `run_summary.json` — machine-readable run report (inputs, parameters, system info)
- `run_summary.txt` — human-readable run report with reproducibility command

## CLI Reference

```
Usage: seq_matcher [OPTIONS] --primer-csv <PRIMER_CSV> --library-csv <LIBRARY_CSV>

Options:
  -p, --primer-csv <PRIMER_CSV>
          Primer CSV file path (columns: id, forward_seq, reverse_seq)
  -l, --library-csv <LIBRARY_CSV>
          Library variant CSV file path
      --library-seq-col <LIBRARY_SEQ_COL>
          Column name in library CSV containing variant sequence
          [default: single_degenerate_library_expanded_reference]
  -s, --seq <SEQ_FILES>
          Sequence files in LABEL:PATH format (e.g. a_11:data/11_seq.txt).
          Repeatable for batch processing.
  -o, --output-dir <OUTPUT_DIR>
          Output directory for result CSV files [default: output]
  -c, --chunk-size <CHUNK_SIZE>
          Sequences per parallel processing chunk [default: 100000]
  -t, --threads <THREADS>
          Number of worker threads (default: all CPU cores)
      --dry-run
          Validate inputs and exit without processing
  -q, --quiet
          Suppress non-error output
  -v, --verbose
          Enable verbose debug output
      --timestamp-output
          Append timestamp to output filenames for reproducibility
  -h, --help
          Print help information
  -V, --version
          Print version information
```

## Algorithm

1. Load primer and library tables; pre-compute reverse complements for all sequences.
2. Build an **Aho-Corasick automaton** encoding all library variants (original + reverse complement) as multi-pattern matching states. Built once at startup and shared via `Arc` across threads.
3. Estimate total lines from file size (no pre-scan), then stream sequences in configurable chunk sizes for Rayon parallel processing.
4. Each sequence is tested against primers with **first-match-wins** semantics.
5. Matched sequences undergo a **single Aho-Corasick scan** to detect all variant hits (deduplicated), replacing per-variant O(V) `contains()` calls.
6. Each thread accumulates counts lock-free; results are merged at chunk boundaries.
7. Real-time progress bar shows completion percentage, throughput, and ETA.
8. Batch output to CSV files with run summary reports.

### Performance

| Optimization | Implementation |
|-------------|----------------|
| Reverse complement | Compile-time `const` 128-byte LUT, single-cycle indexed mapping |
| Variant matching | Aho-Corasick multi-pattern search, O(L+M) instead of O(V*L) |
| Parallel processing | Rayon work-stealing, chunk-level parallelism |
| Cross-thread sharing | `Arc` zero-copy sharing of primers, library, automaton |
| Output I/O | `BufWriter` buffered writes |
| Line counting | File-size-based estimation avoids full pre-scan I/O |
| Variant count output | Pre-built column vectors for O(1) array lookups over O(P*V) HashMap access |
| Progress feedback | `\r` in-place progress bar, no extra I/O overhead |

## Library Usage

SeqMatcher can be used as a Rust library:

```rust
use seq_matcher::{io, matcher};

let primer_data = io::load_primers("primers.csv")?;
let lib = io::load_library("library.csv", "sequence_col")?;
let ac_data = matcher::build_variant_ac(&lib.variants);
```

## Dependencies

| Crate | Purpose |
|-------|---------|
| clap | CLI argument parsing |
| csv | CSV reading and writing |
| rayon | Data parallelism |
| anyhow | Error handling with context |
| aho-corasick | Multi-pattern substring matching (variant search) |
| tracing | Structured logging |
| tracing-subscriber | Log formatting and env-filter support |

### Dev Dependencies

| Crate | Purpose |
|-------|---------|
| criterion | Statistical benchmarking |
| tempfile | Temporary files for integration tests |

## Citation

If you use SeqMatcher in your research, please cite:

```bibtex
@software{seq_matcher,
  author = {Zhao, Jiwen},
  title = {SeqMatcher: High-performance DNA sequence primer matching and variant counting},
  year = {2026},
  version = {2.0.0},
  url = {https://github.com/CropCoder/SeqMatcher}
}
```

A `CITATION.cff` file is included in the repository for GitHub/Zenodo integration.

## License

MIT — see [LICENSE](LICENSE) for details.
