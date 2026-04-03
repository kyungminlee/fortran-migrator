#!/usr/bin/env bash
# 05-compile.sh — Compile migrated BLAS files to verify syntactic correctness.
#
# Compiles each migrated .f/.f90 file individually with gfortran.
# This catches syntax errors introduced by the migration.
#
# Usage: ./recipes/blas/05-compile.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Compiling migrated BLAS (KIND=$TARGET_KIND) ==="

# --------------------------------------------------------------------------
# Find Fortran compiler
# --------------------------------------------------------------------------
FC="${FC:-gfortran}"
if ! command -v "$FC" &>/dev/null; then
    echo "Error: Fortran compiler not found: $FC" >&2
    echo "  Set FC environment variable or install gfortran." >&2
    exit 1
fi

echo "Compiler: $FC ($($FC --version 2>&1 | head -1))"
echo ""

# --------------------------------------------------------------------------
# Compile flags
# --------------------------------------------------------------------------
FFLAGS="-c -O0 -w"
# For fixed-form, gfortran auto-detects from .f extension
# For free-form, gfortran auto-detects from .f90 extension

COMPILE_REPORT="${REPORT_DIR}/compile.txt"
PASS=0
FAIL=0
SKIP=0

{
    echo "# Compilation Report"
    echo "# TARGET_KIND=$TARGET_KIND"
    echo "# Compiler: $FC"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo ""

    if [[ ! -d "$OUTPUT_DIR" ]] || ! ls "$OUTPUT_DIR"/*.f &>/dev/null 2>&1; then
        echo "No output files to compile. Run 02-migrate.sh first."
        exit 1
    fi

    cd "$BUILD_DIR"

    for src in "$OUTPUT_DIR"/*.f "$OUTPUT_DIR"/*.f90; do
        [[ -f "$src" ]] || continue
        base=$(basename "$src")
        obj="${base%.*}.o"

        if $FC $FFLAGS -o "$obj" "$src" 2>"${base}.err"; then
            echo "  PASS: $base"
            ((PASS++))
            rm -f "${base}.err"
        else
            echo "  FAIL: $base"
            sed 's/^/    /' "${base}.err"
            ((FAIL++))
        fi
    done

    echo ""
    echo "# Summary"
    echo "  PASS: $PASS"
    echo "  FAIL: $FAIL"
    echo "  SKIP: $SKIP"

} > "$COMPILE_REPORT" 2>&1

cat "$COMPILE_REPORT"
echo ""
echo "Report: $COMPILE_REPORT"

if [[ "$FAIL" -gt 0 ]]; then
    echo ""
    echo "COMPILATION FAILED: $FAIL file(s) had errors."
    exit 1
else
    echo ""
    echo "ALL FILES COMPILED SUCCESSFULLY."
fi
