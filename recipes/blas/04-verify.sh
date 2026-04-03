#!/usr/bin/env bash
# 04-verify.sh — Verify migrated BLAS files for correctness.
#
# Checks:
#   1. Output files exist (at least as many as expected)
#   2. No residual source-precision types in code lines
#   3. No residual D-exponent literals in code lines
#   4. Column width (fixed-form lines ≤ 72 chars)
#
# Usage: ./recipes/blas/04-verify.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Verifying migrated BLAS (KIND=$TARGET_KIND) ==="

VERIFY_REPORT="${REPORT_DIR}/verify.txt"
ERRORS=0
WARNINGS=0

{
    echo "# Verification Report"
    echo "# TARGET_KIND=$TARGET_KIND"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo ""

    # ------------------------------------------------------------------
    # Check 1: Output files exist
    # ------------------------------------------------------------------
    echo "## Check 1: Output file count"
    echo ""

    if [[ -d "$OUTPUT_DIR" ]]; then
        out_count=$(find "$OUTPUT_DIR" -maxdepth 1 \( -name '*.f' -o -name '*.f90' \) | wc -l | tr -d ' ')
        echo "  Output files: $out_count"
        if [[ "$out_count" -lt 70 ]]; then
            echo "  ERROR: Expected at least 70 output files, got $out_count"
            ((ERRORS++))
        else
            echo "  OK: $out_count files (≥ 70 expected)"
        fi
        # Check utility files exist
        for f in lsame.f xerbla.f xerbla_array.f; do
            if [[ -f "${OUTPUT_DIR}/${f}" ]]; then
                echo "  OK: $f"
            else
                echo "  MISSING: $f"
                ((ERRORS++))
            fi
        done
    else
        echo "  ERROR: Output directory does not exist: $OUTPUT_DIR"
        ((ERRORS++))
    fi

    echo ""

    # ------------------------------------------------------------------
    # Check 2: No residual precision types in code (non-comment) lines
    # ------------------------------------------------------------------
    echo "## Check 2: Residual precision types"
    echo ""

    if [[ -d "$OUTPUT_DIR" ]] && ls "$OUTPUT_DIR"/*.f &>/dev/null; then
        for f in "$OUTPUT_DIR"/*.f; do
            base=$(basename "$f")
            # Skip precision-independent files
            [[ "$base" == "lsame.f" || "$base" == "xerbla.f" || "$base" == "xerbla_array.f" ]] && continue
            # Skip i-prefix files (integer return, may legitimately have typed args)
            [[ "$base" == i*.f ]] && continue

            # Look for DOUBLE PRECISION, COMPLEX*16 in non-comment code lines
            while IFS= read -r line; do
                echo "  RESIDUAL TYPE in $base: $line"
                ((WARNINGS++))
            done < <(grep -n '^[^Cc*!]' "$f" 2>/dev/null \
                | grep -iE 'DOUBLE PRECISION|COMPLEX\*16|COMPLEX\*8|DOUBLE COMPLEX' || true)
        done

        # Free-form: check for old wp definitions
        for f in "$OUTPUT_DIR"/*.f90; do
            [[ -f "$f" ]] || continue
            while IFS= read -r line; do
                echo "  RESIDUAL WP in $(basename "$f"): $line"
                ((WARNINGS++))
            done < <(grep -n 'kind(1\.[de]0)' "$f" 2>/dev/null || true)
        done
    else
        echo "  SKIP: No output files to check"
    fi

    [[ "$WARNINGS" -eq 0 ]] && echo "  OK: No residual types found"
    echo ""

    # ------------------------------------------------------------------
    # Check 3: No residual D-exponent literals in code lines
    # ------------------------------------------------------------------
    echo "## Check 3: Residual D-exponent literals"
    echo ""

    lit_warnings=0
    if [[ -d "$OUTPUT_DIR" ]] && ls "$OUTPUT_DIR"/*.f &>/dev/null; then
        for f in "$OUTPUT_DIR"/*.f; do
            base=$(basename "$f")
            [[ "$base" == "lsame.f" || "$base" == "xerbla.f" || "$base" == "xerbla_array.f" ]] && continue
            while IFS= read -r line; do
                echo "  RESIDUAL D-LITERAL in $base: $line"
                ((lit_warnings++))
            done < <(grep -n '^[^Cc*!]' "$f" 2>/dev/null \
                | grep -iE '[0-9]\.[0-9]*[Dd][+-]?[0-9]' || true)
        done
    else
        echo "  SKIP: No output files to check"
    fi

    [[ "$lit_warnings" -eq 0 ]] && echo "  OK: No residual D-exponent literals"
    echo ""

    # ------------------------------------------------------------------
    # Check 4: Column width for fixed-form files
    # ------------------------------------------------------------------
    echo "## Check 4: Fixed-form column width (max 72)"
    echo ""

    col_warnings=0
    if [[ -d "$OUTPUT_DIR" ]] && ls "$OUTPUT_DIR"/*.f &>/dev/null; then
        for f in "$OUTPUT_DIR"/*.f; do
            # Only check non-comment lines (column 72 rule doesn't apply to comments)
            while IFS= read -r line; do
                echo "  $line"
                ((col_warnings++))
            done < <(awk '/^[^Cc*!]/ && length > 72 { printf "OVERFLOW in %s line %d: %d chars\n", FILENAME, NR, length }' "$f")
        done
    else
        echo "  SKIP: No output files to check"
    fi

    [[ "$col_warnings" -eq 0 ]] && echo "  OK: All lines within 72 columns"
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
