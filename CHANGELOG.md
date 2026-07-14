# Changelog

All notable changes to SeqMatcher will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] - 2026-07-14

### Added
- `build-release.sh` cross-platform release builder (macOS/Linux/Windows targets)
- `--no-report` flag to suppress run summary report generation
- Pre-built release binaries for macOS (arm64/x86_64) and Linux (amd64/arm64 musl)

### Changed
- Simplified `progress` module: removed unused `count_lines` public API in favor of `estimate_total_lines`
- Renamed internal `count_lines_in_buf` to `count_lines_exact` for clarity
- Updated `rustfmt.toml` to use only stable-channel formatting options
- `ThreadResult` now derives `Default` instead of manual implementation

### Fixed
- Clippy `too_many_arguments` warning on `process_sequences` (annotated allowed)
- Clippy `derivable_impls` warning on `ThreadResult::Default`
- Unused `create_temp_seq` helper in integration tests (annotated reserved)

## [2.0.0] - 2026-07-14

### Added
- Comprehensive IUPAC degenerate base support (S, W, B, D, H, V) in reverse complement table
- GitHub Actions CI pipeline with format/lint/test/build matrix
- Dependabot configuration for automated dependency updates
- Makefile with standard targets (build, test, lint, fmt, bench, clean)
- `--dry-run` mode for input validation without processing
- `--quiet` and `--verbose` output control flags
- `--timestamp-output` flag for time-stamped output filenames
- Input data validation with health checks
- Run summary report (`run_summary.json` and `run_summary.txt`)
- Reproducibility features: full CLI arg recording in summary reports
- `CITATION.cff` for academic citation
- Benchmark suite using Criterion (AC matching, reverse complement)
- Integration tests with end-to-end pipeline coverage
- Public library API via `src/lib.rs` for downstream crate consumption
- Comprehensive rustdoc documentation for all public APIs

### Changed
- Refactored `main.rs` (436 lines) into focused modules: `cli`, `progress`, `matcher`, `pipeline`
- Replaced manual `csv_to_line` serialization with `csv::Writer` for correctness
- Optimized `write_variant_counts` O(P*V) HashMap lookups to O(1) array indexing
- Replaced `eprintln!` diagnostics with `tracing` structured logging
- Unified error messages to English (CLI help text remains bilingual)
- Eliminated double I/O from `count_lines` pre-scan; uses file-size estimation

### Fixed
- **Critical**: IUPAC degenerate bases S, W, B, D, H, V now correctly reverse-complemented
  (previously mapped to themselves, producing incorrect results for degenerate sequences)
- Version inconsistency between `Cargo.toml` (1.4.0) and `Cargo.lock` (1.3.0)

### Removed
- Unused `serde` dependency from `Cargo.toml`
- Manual `csv_to_line` function in favor of `csv::Writer`

## [1.4.0] - 2025-06-25

### Added
- Streaming I/O with reusable hit buffer for reduced allocations
- In-place ASCII uppercase conversion

## [1.3.0] - 2025-06-25

### Added
- Aho-Corasick multi-pattern matching for variant search (replaces O(V) naive scan)

## [1.2.0] - 2025-06-24

### Added
- Progress bar with percentage, throughput rate, and ETA display

## [1.1.0] - 2025-06-24

### Added
- `build-linux.sh` cross-compilation script for Linux musl target

## [1.0.0] - 2025-06-24

### Added
- Initial release: high-performance DNA sequence primer matching CLI
- CSV-based primer and library variant loading
- Rayon-powered parallel chunk processing
- Reverse complement computation with compile-time LUT
- First-match-wins primer matching with prefix/suffix anchoring
