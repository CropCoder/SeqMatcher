#!/bin/bash
# =============================================================================
# SeqMatcher cross-platform release builder
# Builds release binaries for Linux, macOS, and Windows from macOS host.
# Requires: rustup with all targets installed
# Author: Jiwen Zhao (https://github.com/CropCoder)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOOLCHAIN_DIR="$(rustup show home)/toolchains/$(rustup show active-toolchain | awk '{print $1}')"
RELEASE_DIR="$SCRIPT_DIR/releases"
VERSION=$(grep '^version' "$SCRIPT_DIR/Cargo.toml" | head -1 | sed 's/.*"\(.*\)"/\1/')
TIMESTAMP=$(date -u +"%Y%m%d-%H%M%S")

# -- color support --
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # no color

# ---------------------------------------------------------------------------
# Platform definitions
#   target-triple | binary-name-on-disk | display-name | os-tag
# ---------------------------------------------------------------------------
declare -A TARGETS
TARGETS=(
    ["x86_64-unknown-linux-musl"]="seq_matcher|Linux x86_64 (musl static)|linux-amd64"
    ["aarch64-unknown-linux-musl"]="seq_matcher|Linux aarch64 (musl static)|linux-arm64"
    ["x86_64-apple-darwin"]="seq_matcher|macOS x86_64|macos-amd64"
    ["aarch64-apple-darwin"]="seq_matcher|macOS ARM64 (Apple Silicon)|macos-arm64"
    ["x86_64-pc-windows-gnu"]="seq_matcher.exe|Windows x86_64|windows-amd64"
)

echo -e "${GREEN}=== SeqMatcher v${VERSION} Cross-Platform Release Builder ===${NC}"
echo "Timestamp: $TIMESTAMP"
echo "Rust toolchain: $(rustup show active-toolchain)"
echo ""

# -- Ensure all targets are installed --
echo -e "${YELLOW}[1/4] Checking target toolchains...${NC}"
for target in "${!TARGETS[@]}"; do
    if [ ! -d "$TOOLCHAIN_DIR/lib/rustlib/$target" ]; then
        echo "  Installing $target ..."
        rustup target add "$target"
    else
        echo "  $target ... OK"
    fi
done
echo ""

# -- Build all targets --
echo -e "${YELLOW}[2/4] Building release binaries...${NC}"
mkdir -p "$RELEASE_DIR"

for target in "${!TARGETS[@]}"; do
    IFS='|' read -r binary_name display_name os_tag <<< "${TARGETS[$target]}"
    echo "  Building for $display_name ($target) ..."

    # Windows GNU target needs a cross-linker on macOS
    if [ "$target" = "x86_64-pc-windows-gnu" ] && ! command -v x86_64-w64-mingw32-gcc &>/dev/null; then
        echo -e "    ${RED}ERROR: x86_64-w64-mingw32-gcc not found. Install with: brew install mingw-w64${NC}"
        continue
    fi

    # Linux musl targets need musl-cross on macOS
    if [[ "$target" == *linux*musl* ]]; then
        case "$target" in
            x86_64-unknown-linux-musl)
                if ! command -v x86_64-linux-musl-gcc &>/dev/null; then
                    echo -e "    ${RED}ERROR: x86_64-linux-musl-gcc not found. Install with: brew install musl-cross${NC}"
                    continue
                fi
                ;;
            aarch64-unknown-linux-musl)
                if ! command -v aarch64-linux-musl-gcc &>/dev/null; then
                    echo -e "    ${YELLOW}    WARNING: aarch64-linux-musl-gcc not found. Skipping.${NC}"
                    echo -e "    ${YELLOW}    Install with: brew install aarch64-unknown-linux-musl-binutils${NC}"
                    continue
                fi
                ;;
        esac
    fi

    # Use env -i for musl targets to avoid conda/brew compiler pollution
    if [[ "$target" == *linux*musl* ]]; then
        env -i \
            HOME="$HOME" \
            PATH="$TOOLCHAIN_DIR/bin:/opt/homebrew/bin:/usr/bin:/bin" \
            RUSTUP_HOME="$(rustup show home)" \
            CARGO_HOME="${CARGO_HOME:-$HOME/.cargo}" \
            cargo build --release --target "$target" --manifest-path "$SCRIPT_DIR/Cargo.toml"
    else
        cargo build --release --target "$target"
    fi

    # Copy binary to releases dir with descriptive name
    src="$SCRIPT_DIR/target/$target/release/$binary_name"
    dst="$RELEASE_DIR/seq_matcher-v${VERSION}-${os_tag}"
    if [ -f "$src" ]; then
        cp "$src" "$dst"
        chmod +x "$dst" 2>/dev/null || true
        size=$(ls -lh "$dst" | awk '{print $5}')
        echo -e "  ${GREEN}  -> $dst ($size)${NC}"
    else
        echo -e "  ${RED}  ERROR: binary not found at $src${NC}"
    fi
done
echo ""

# -- Create checksums --
echo -e "${YELLOW}[3/4] Generating SHA256 checksums...${NC}"
CHECKSUM_FILE="$RELEASE_DIR/seq_matcher-v${VERSION}-checksums-sha256.txt"
> "$CHECKSUM_FILE"
for file in "$RELEASE_DIR"/seq_matcher-v${VERSION}-*; do
    if [ -f "$file" ] && [[ "$file" != *.txt ]]; then
        basename=$(basename "$file")
        shasum -a 256 "$file" | sed "s|$file|$basename|" >> "$CHECKSUM_FILE"
    fi
done
cat "$CHECKSUM_FILE"
echo ""

# -- Generate release manifest --
echo -e "${YELLOW}[4/4] Generating release manifest...${NC}"
MANIFEST="$RELEASE_DIR/seq_matcher-v${VERSION}-manifest.txt"
{
    echo "SeqMatcher v${VERSION} Release Manifest"
    echo "======================================="
    echo "Build date: $(date -u)"
    echo "Git commit: $(git -C "$SCRIPT_DIR" rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
    echo "Rust toolchain: $(rustup show active-toolchain)"
    echo ""
    echo "Binaries:"
    for file in "$RELEASE_DIR"/seq_matcher-v${VERSION}-*; do
        if [ -f "$file" ] && [[ "$file" != *.txt ]]; then
            basename=$(basename "$file")
            echo "  $basename"
        fi
    done
    echo ""
    echo "Installation:"
    echo "  chmod +x <binary> && ./<binary> --help"
    echo ""
    echo "Verification:"
    echo "  shasum -a 256 -c seq_matcher-v${VERSION}-checksums-sha256.txt"
} > "$MANIFEST"
cat "$MANIFEST"

echo ""
echo -e "${GREEN}=== Release complete! ===${NC}"
echo -e "Output directory: ${BLUE}$RELEASE_DIR${NC}"
echo ""
ls -lh "$RELEASE_DIR"
