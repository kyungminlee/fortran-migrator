#!/usr/bin/env bash
# 06-build-library.sh — Build migrated BLAS into a static library.
#
# Compiles all migrated source files and archives them into a static
# library (e.g., libblas_q.a for KIND=16).
#
# Usage: ./recipes/blas/06-build-library.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Building BLAS library (KIND=$TARGET_KIND) ==="

# --------------------------------------------------------------------------
# Determine library name
# --------------------------------------------------------------------------
if [[ "$TARGET_KIND" == "10" ]]; then
    LIB_NAME="libblas_e.a"
else
    LIB_NAME="libblas_q.a"
fi

LIB_PATH="${BUILD_DIR}/${LIB_NAME}"

# --------------------------------------------------------------------------
# Find tools
# --------------------------------------------------------------------------
FC="${FC:-gfortran}"
AR="${AR:-ar}"
RANLIB="${RANLIB:-ranlib}"

if ! command -v "$FC" &>/dev/null; then
    echo "Error: Fortran compiler not found: $FC" >&2
    exit 1
fi

echo "Compiler: $FC"
echo "Archiver: $AR"
echo "Library:  $LIB_PATH"
echo ""

# --------------------------------------------------------------------------
# Compile flags
# --------------------------------------------------------------------------
FFLAGS="${FFLAGS:--O2 -w}"

# --------------------------------------------------------------------------
# Compile
# --------------------------------------------------------------------------
if [[ ! -d "$OUTPUT_DIR" ]] || ! ls "$OUTPUT_DIR"/*.f &>/dev/null 2>&1; then
    echo "Error: No migrated files found in $OUTPUT_DIR" >&2
    echo "  Run 02-migrate.sh first." >&2
    exit 1
fi

cd "$BUILD_DIR"
rm -f *.o "$LIB_NAME"

SOURCES=()
for src in "$OUTPUT_DIR"/*.f "$OUTPUT_DIR"/*.f90; do
    [[ -f "$src" ]] || continue
    SOURCES+=("$src")
done

echo "Compiling ${#SOURCES[@]} files..."

FAIL=0
OBJECTS=()
for src in "${SOURCES[@]}"; do
    base=$(basename "$src")
    obj="${base%.*}.o"
    if $FC $FFLAGS -c -o "$obj" "$src" 2>/dev/null; then
        OBJECTS+=("$obj")
    else
        echo "  FAIL: $base"
        ((FAIL++))
    fi
done

if [[ "$FAIL" -gt 0 ]]; then
    echo ""
    echo "Error: $FAIL file(s) failed to compile." >&2
    echo "  Run 05-compile.sh for detailed error messages." >&2
    exit 1
fi

# --------------------------------------------------------------------------
# Archive
# --------------------------------------------------------------------------
echo "Archiving ${#OBJECTS[@]} objects into $LIB_NAME..."
$AR rcs "$LIB_NAME" "${OBJECTS[@]}"
$RANLIB "$LIB_NAME"

# --------------------------------------------------------------------------
# Symbol verification
# --------------------------------------------------------------------------
SYMBOL_REPORT="${REPORT_DIR}/library-symbols.txt"

echo ""
echo "=== Library Symbol Report ==="

nm "$LIB_PATH" | grep ' [TtWw] ' | awk '{print $3}' | sort > "$SYMBOL_REPORT"

total=$(wc -l < "$SYMBOL_REPORT" | tr -d ' ')
echo "Total exported symbols: $total"
echo "Symbol list: $SYMBOL_REPORT"

# Show a sample
echo ""
echo "Sample symbols:"
head -20 "$SYMBOL_REPORT" | sed 's/^/  /'
if [[ "$total" -gt 20 ]]; then
    echo "  ... ($((total - 20)) more)"
fi

echo ""
echo "Library built: $LIB_PATH ($(du -h "$LIB_PATH" | cut -f1))"
