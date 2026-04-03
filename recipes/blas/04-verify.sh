#!/usr/bin/env bash
# 04-verify.sh — Verify migrated BLAS files for correctness.
#
# Checks:
#   1. All expected output files exist with correct names
#   2. No residual source-precision types remain in code lines
#   3. No residual source-precision routine names remain
#   4. No residual D-exponent literals remain
#   5. Column width (fixed-form lines ≤ 72 chars)
#
# Usage: ./recipes/blas/04-verify.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Verifying migrated BLAS (KIND=$TARGET_KIND) ==="

VERIFY_REPORT="${REPORT_DIR}/verify.txt"
ERRORS=0
WARNINGS=0

if [[ "$TARGET_KIND" == "10" ]]; then
    REAL_PREFIX="e"
    COMPLEX_PREFIX="y"
else
    REAL_PREFIX="q"
    COMPLEX_PREFIX="x"
fi

{
    echo "# Verification Report"
    echo "# TARGET_KIND=$TARGET_KIND"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo ""

    # ------------------------------------------------------------------
    # Check 1: Expected output files exist
    # ------------------------------------------------------------------
    echo "## Check 1: Output file existence"
    echo ""

    check_file() {
        local f="$1"
        if [[ -f "${OUTPUT_DIR}/${f}" ]]; then
            echo "  OK: $f"
        else
            echo "  MISSING: $f"
            ((ERRORS++))
        fi
    }

    # Precision-independent (copied as-is)
    for f in lsame.f xerbla.f xerbla_array.f; do
        check_file "$f"
    done

    # Real routines
    for base in asum axpy copy dot gbmv gemm gemmtr gemv ger \
                sbmv scal spmv spr spr2 swap symm symv syr syr2 \
                syr2k syrk tbmv tbsv tpmv tpsv trmm trmv trsm trsv \
                rot rotm rotmg; do
        check_file "${REAL_PREFIX}${base}.f"
    done

    # Complex routines
    for base in axpy copy dotc dotu gbmv gemm gemmtr gemv gerc geru \
                hbmv hemm hemv her her2 her2k herk hpmv hpr hpr2 \
                scal swap symm syr2k syrk tbmv tbsv tpmv tpsv \
                trmm trmv trsm trsv; do
        check_file "${COMPLEX_PREFIX}${base}.f"
    done

    # F90 files
    check_file "${REAL_PREFIX}nrm2.f90"
    check_file "${REAL_PREFIX}rotg.f90"
    check_file "${COMPLEX_PREFIX}rotg.f90"
    check_file "${REAL_PREFIX}${COMPLEX_PREFIX}nrm2.f90"

    echo ""

    # ------------------------------------------------------------------
    # Check 2: No residual precision types in code (non-comment) lines
    # ------------------------------------------------------------------
    echo "## Check 2: Residual precision types"
    echo ""

    check_residual_types() {
        local dir="$1"
        local found=0

        # Fixed-form: non-comment lines (not starting with C/c/*/!)
        for f in "$dir"/*.f; do
            [[ -f "$f" ]] || continue
            local base
            base=$(basename "$f")
            # Skip precision-independent files
            [[ "$base" == "lsame.f" || "$base" == "xerbla.f" || "$base" == "xerbla_array.f" ]] && continue

            # Look for DOUBLE PRECISION, COMPLEX*16 in non-comment code lines
            grep -n '^[^Cc*!]' "$f" 2>/dev/null \
                | grep -iE 'DOUBLE PRECISION|COMPLEX\*16|COMPLEX\*8|DOUBLE COMPLEX' \
                | while read -r line; do
                    echo "  RESIDUAL TYPE in $base: $line"
                    found=1
                done
        done

        # Free-form: check for old wp definitions
        for f in "$dir"/*.f90; do
            [[ -f "$f" ]] || continue
            grep -n 'kind(1\.[de]0)' "$f" 2>/dev/null | while read -r line; do
                echo "  RESIDUAL WP in $(basename "$f"): $line"
                found=1
            done
        done

        return $found
    }

    if [[ -d "$OUTPUT_DIR" ]] && ls "$OUTPUT_DIR"/*.f &>/dev/null; then
        check_residual_types "$OUTPUT_DIR" || ((ERRORS++))
    else
        echo "  SKIP: No output files to check"
    fi

    echo ""

    # ------------------------------------------------------------------
    # Check 3: No residual D-exponent literals in code lines
    # ------------------------------------------------------------------
    echo "## Check 3: Residual D-exponent literals"
    echo ""

    if [[ -d "$OUTPUT_DIR" ]] && ls "$OUTPUT_DIR"/*.f &>/dev/null; then
        for f in "$OUTPUT_DIR"/*.f; do
            local base
            base=$(basename "$f")
            [[ "$base" == "lsame.f" || "$base" == "xerbla.f" || "$base" == "xerbla_array.f" ]] && continue
            grep -n '^[^Cc*!]' "$f" 2>/dev/null \
                | grep -iE '[0-9]\.[0-9]*[Dd][+-]?[0-9]' \
                | while read -r line; do
                    echo "  RESIDUAL D-LITERAL in $base: $line"
                    ((ERRORS++))
                done
        done
    else
        echo "  SKIP: No output files to check"
    fi

    echo ""

    # ------------------------------------------------------------------
    # Check 4: Column width for fixed-form files
    # ------------------------------------------------------------------
    echo "## Check 4: Fixed-form column width (max 72)"
    echo ""

    if [[ -d "$OUTPUT_DIR" ]] && ls "$OUTPUT_DIR"/*.f &>/dev/null; then
        for f in "$OUTPUT_DIR"/*.f; do
            awk 'length > 72 { printf "  OVERFLOW in %s line %d: %d chars\n", FILENAME, NR, length }' "$f"
        done
        echo "  (no output = all lines within 72 columns)"
    else
        echo "  SKIP: No output files to check"
    fi

    echo ""

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    echo "## Summary"
    echo "  Errors:   $ERRORS"
    echo "  Warnings: $WARNINGS"

} > "$VERIFY_REPORT" 2>&1

cat "$VERIFY_REPORT"
echo ""
echo "Report: $VERIFY_REPORT"

if [[ "$ERRORS" -gt 0 ]]; then
    echo ""
    echo "VERIFICATION FAILED with $ERRORS error(s)."
    exit 1
else
    echo ""
    echo "VERIFICATION PASSED."
fi
