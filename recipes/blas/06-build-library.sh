#!/usr/bin/env bash
# 06-build-library.sh — Build migrated BLAS into static libraries.
#
# Produces two libraries:
#   libblas_common.a  — type-independent routines (LSAME, XERBLA, ...)
#   libqblas.a        — precision-specific routines (Q/X/I-prefix)
#
# Usage: ./recipes/blas/06-build-library.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Building BLAS libraries (KIND=$TARGET_KIND) ==="

# --------------------------------------------------------------------------
# Determine library names
# --------------------------------------------------------------------------
if [[ "$TARGET_KIND" == "10" ]]; then
    LIB_NAME="libeblas.a"
else
    LIB_NAME="libqblas.a"
fi
COMMON_LIB_NAME="libblas_common.a"

# Type-independent source files (no precision prefix)
COMMON_FILES="lsame.f xerbla.f xerbla_array.f"

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
echo ""

# --------------------------------------------------------------------------
# Compile flags
# --------------------------------------------------------------------------
FFLAGS="${FFLAGS:--O2 -w}"

# --------------------------------------------------------------------------
# Check input
# --------------------------------------------------------------------------
if [[ ! -d "$OUTPUT_DIR" ]] || ! ls "$OUTPUT_DIR"/*.f &>/dev/null 2>&1; then
    echo "Error: No migrated files found in $OUTPUT_DIR" >&2
    echo "  Run 02-migrate.sh first." >&2
    exit 1
fi

cd "$BUILD_DIR"
rm -f *.o *.a

# --------------------------------------------------------------------------
# Classify source files into common vs precision-specific
# --------------------------------------------------------------------------
COMMON_SOURCES=()
PRECISION_SOURCES=()

for src in "$OUTPUT_DIR"/*.f "$OUTPUT_DIR"/*.f90; do
    [[ -f "$src" ]] || continue
    base=$(basename "$src")
    is_common=false
    for cf in $COMMON_FILES; do
        if [[ "$base" == "$cf" ]]; then
            is_common=true
            break
        fi
    done
    if $is_common; then
        COMMON_SOURCES+=("$src")
    else
        PRECISION_SOURCES+=("$src")
    fi
done

echo "Common files:    ${#COMMON_SOURCES[@]}"
echo "Precision files: ${#PRECISION_SOURCES[@]}"
echo ""

# --------------------------------------------------------------------------
# Compile and archive common library
# --------------------------------------------------------------------------
echo "--- Building $COMMON_LIB_NAME ---"
COMMON_OBJECTS=()
FAIL=0
for src in "${COMMON_SOURCES[@]}"; do
    base=$(basename "$src")
    obj="${base%.*}.o"
    if $FC $FFLAGS -c -o "$obj" "$src" 2>/dev/null; then
        COMMON_OBJECTS+=("$obj")
    else
        echo "  FAIL: $base"
        FAIL=$((FAIL+1))
    fi
done

if [[ ${#COMMON_OBJECTS[@]} -gt 0 ]]; then
    $AR rcs "$COMMON_LIB_NAME" "${COMMON_OBJECTS[@]}"
    $RANLIB "$COMMON_LIB_NAME"
    echo "  Built: $COMMON_LIB_NAME ($(du -h "$COMMON_LIB_NAME" | cut -f1))"
    nm "$COMMON_LIB_NAME" | grep ' T ' | awk '{print $3}' | sed 's/^_//;s/_$//' | tr 'a-z' 'A-Z' | grep -v LTMP | sort -u | sed 's/^/    /'
fi
echo ""

# --------------------------------------------------------------------------
# Compile and archive precision-specific library
# --------------------------------------------------------------------------
echo "--- Building $LIB_NAME ---"
PRECISION_OBJECTS=()
for src in "${PRECISION_SOURCES[@]}"; do
    base=$(basename "$src")
    obj="prec_${base%.*}.o"  # prefix to avoid name collision with common
    if $FC $FFLAGS -c -o "$obj" "$src" 2>/dev/null; then
        PRECISION_OBJECTS+=("$obj")
    else
        echo "  FAIL: $base"
        FAIL=$((FAIL+1))
    fi
done

if [[ "$FAIL" -gt 0 ]]; then
    echo ""
    echo "Error: $FAIL file(s) failed to compile." >&2
    echo "  Run 05-compile.sh for detailed error messages." >&2
    exit 1
fi

$AR rcs "$LIB_NAME" "${PRECISION_OBJECTS[@]}"
$RANLIB "$LIB_NAME"
echo "  Built: $LIB_NAME ($(du -h "$LIB_NAME" | cut -f1))"

# --------------------------------------------------------------------------
# Symbol report
# --------------------------------------------------------------------------
SYMBOL_REPORT="${REPORT_DIR}/library-symbols.txt"

{
    echo "# Library Symbol Report"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo ""
    echo "## $COMMON_LIB_NAME"
    nm "$COMMON_LIB_NAME" | grep ' T ' | awk '{print $3}' | sed 's/^_//;s/_$//' | tr 'a-z' 'A-Z' | grep -v LTMP | sort -u
    echo ""
    echo "## $LIB_NAME"
    nm "$LIB_NAME" | grep ' T ' | awk '{print $3}' | sed 's/^_//;s/_$//' | tr 'a-z' 'A-Z' | grep -v LTMP | sort -u
} > "$SYMBOL_REPORT"

q_count=$(nm "$LIB_NAME" | grep ' T ' | awk '{print $3}' | sed 's/^_//;s/_$//' | tr 'a-z' 'A-Z' | grep '^Q' | sort -u | wc -l | tr -d ' ')
x_count=$(nm "$LIB_NAME" | grep ' T ' | awk '{print $3}' | sed 's/^_//;s/_$//' | tr 'a-z' 'A-Z' | grep '^X' | grep -v XERBLA | sort -u | wc -l | tr -d ' ')
i_count=$(nm "$LIB_NAME" | grep ' T ' | awk '{print $3}' | sed 's/^_//;s/_$//' | tr 'a-z' 'A-Z' | grep '^I' | sort -u | wc -l | tr -d ' ')

echo ""
echo "=== Summary ==="
echo "  $COMMON_LIB_NAME: type-independent routines"
echo "  $LIB_NAME:        $q_count Q-prefix + $x_count X-prefix + $i_count I-prefix"
echo ""
echo "Symbol report: $SYMBOL_REPORT"
