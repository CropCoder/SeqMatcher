#!/bin/bash
# Build Linux x86_64 static binary (musl)
# Works even when conda overrides rustup's rustc

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOOLCHAIN_DIR="$(rustup show home)/toolchains/$(rustup show active-toolchain | awk '{print $1}')"

echo "==> Building SeqMatcher for Linux x86_64 (musl static)..."

if [ ! -f "$TOOLCHAIN_DIR/bin/cargo" ]; then
    echo "ERROR: rustup toolchain not found at $TOOLCHAIN_DIR"
    exit 1
fi

if [ ! -d "$TOOLCHAIN_DIR/lib/rustlib/x86_64-unknown-linux-musl" ]; then
    echo "==> Installing x86_64-unknown-linux-musl target..."
    rustup target add x86_64-unknown-linux-musl
fi

if ! command -v x86_64-linux-musl-gcc &>/dev/null; then
    echo "==> Installing musl-cross linker..."
    brew install musl-cross
fi

env -i \
    HOME="$HOME" \
    PATH="$TOOLCHAIN_DIR/bin:/opt/homebrew/bin:/usr/bin:/bin" \
    RUSTUP_HOME="$(rustup show home)" \
    CARGO_HOME="${CARGO_HOME:-$HOME/.cargo}" \
    cargo build --release --target x86_64-unknown-linux-musl --manifest-path "$SCRIPT_DIR/Cargo.toml"

echo ""
echo "==> Done: $SCRIPT_DIR/target/x86_64-unknown-linux-musl/release/seq_matcher"
file "$SCRIPT_DIR/target/x86_64-unknown-linux-musl/release/seq_matcher"
