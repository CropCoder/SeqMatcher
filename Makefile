.PHONY: build release test lint fmt clean check bench

# Default target: debug build
build:
	cargo build

# Optimized release build
release:
	cargo build --release

# Run all tests
test:
	cargo test --release

# Run clippy lint with all warnings as errors
lint:
	cargo clippy -- -D warnings

# Format all source code
fmt:
	cargo fmt -- --check

# Auto-fix formatting
fmt-fix:
	cargo fmt

# Check compilation without producing binary (faster than build)
check:
	cargo check

# Run benchmarks
bench:
	cargo bench

# Clean build artifacts
clean:
	cargo clean

# Full CI pipeline (fmt + lint + test + build)
ci: fmt lint test release
	@echo "All checks passed."
