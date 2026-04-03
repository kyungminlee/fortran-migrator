#!/usr/bin/env bash
# 03-verify.sh — Verify migrated BLACS files for correctness.
#
# Checks:
#   1. All expected cloned files exist
#   2. No residual 'double' type in q-variant files (should be QREAL)
#   3. No residual 'DCOMPLEX' in x-variant files (should be QCOMPLEX)
#   4. No residual MPI_DOUBLE in extended-precision files
#   5. Function names are correctly prefixed
#
# Usage: ./recipes/blacs/03-verify.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Verifying migrated BLACS (KIND=$TARGET_KIND) ==="

RP="$REAL_PREFIX"
CP="$COMPLEX_PREFIX"
ERRORS=0

VERIFY_REPORT="${REPORT_DIR}/verify.txt"

{
    echo "# BLACS Verification Report"
    echo "# TARGET_KIND=$TARGET_KIND"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo ""

    # ------------------------------------------------------------------
    # Check 1: Expected files exist
    # ------------------------------------------------------------------
    echo "## Check 1: Expected output files"

    check_file() {
        if [[ -f "${OUTPUT_DIR}/$1" ]]; then
            echo "  OK: $1"
        else
            echo "  MISSING: $1"
            ((ERRORS++))
        fi
    }

    # Real extended user-facing
    for op in gesd2d_ gerv2d_ gebs2d_ gebr2d_ trsd2d_ trrv2d_ trbs2d_ trbr2d_ gsum2d_ gamx2d_ gamn2d_; do
        check_file "${RP}${op}.c"
    done

    # Complex extended user-facing
    for op in gesd2d_ gerv2d_ gebs2d_ gebr2d_ trsd2d_ trrv2d_ trbs2d_ trbr2d_ gsum2d_ gamx2d_ gamn2d_; do
        check_file "${CP}${op}.c"
    done

    # Real extended helpers
    for f in BI_${RP}vvsum.c BI_${RP}vvamn.c BI_${RP}vvamx.c BI_${RP}mvcopy.c BI_${RP}vmcopy.c; do
        check_file "$f"
    done

    # Complex extended helpers
    for f in BI_${CP}vvsum.c BI_${CP}vvamn.c BI_${CP}vvamx.c; do
        check_file "$f"
    done

    echo ""

    # ------------------------------------------------------------------
    # Check 2: No residual 'double' in q/e-variant files
    # ------------------------------------------------------------------
    echo "## Check 2: Residual 'double' in real-extended files"

    for f in "$OUTPUT_DIR"/${RP}*.c "$OUTPUT_DIR"/BI_${RP}*.c; do
        [[ -f "$f" ]] || continue
        base=$(basename "$f")
        residual=$(grep -c '\bdouble\b' "$f" 2>/dev/null || true)
        if [[ "$residual" -gt 0 ]]; then
            echo "  WARNING: $base has $residual 'double' occurrences"
            grep -n '\bdouble\b' "$f" | head -3 | sed 's/^/    /'
        fi
    done

    echo ""

    # ------------------------------------------------------------------
    # Check 3: No residual DCOMPLEX in x/y-variant files
    # ------------------------------------------------------------------
    echo "## Check 3: Residual 'DCOMPLEX' in complex-extended files"

    for f in "$OUTPUT_DIR"/${CP}*.c "$OUTPUT_DIR"/BI_${CP}*.c; do
        [[ -f "$f" ]] || continue
        base=$(basename "$f")
        residual=$(grep -c '\bDCOMPLEX\b' "$f" 2>/dev/null || true)
        if [[ "$residual" -gt 0 ]]; then
            echo "  WARNING: $base has $residual 'DCOMPLEX' occurrences"
            grep -n '\bDCOMPLEX\b' "$f" | head -3 | sed 's/^/    /'
        fi
    done

    echo ""

    # ------------------------------------------------------------------
    # Check 4: No residual MPI_DOUBLE in extended files
    # ------------------------------------------------------------------
    echo "## Check 4: Residual MPI_DOUBLE in extended files"

    for f in "$OUTPUT_DIR"/${RP}*.c "$OUTPUT_DIR"/${CP}*.c \
             "$OUTPUT_DIR"/BI_${RP}*.c "$OUTPUT_DIR"/BI_${CP}*.c; do
        [[ -f "$f" ]] || continue
        base=$(basename "$f")
        residual=$(grep -c 'MPI_DOUBLE\b' "$f" 2>/dev/null || true)
        if [[ "$residual" -gt 0 ]]; then
            echo "  WARNING: $base has $residual 'MPI_DOUBLE' occurrences"
            grep -n 'MPI_DOUBLE\b' "$f" | head -3 | sed 's/^/    /'
        fi
    done

    echo ""

    # ------------------------------------------------------------------
    # Check 5: Bdef.h has extended types
    # ------------------------------------------------------------------
    echo "## Check 5: Bdef.h extended-precision definitions"

    if grep -q 'QREAL' "$OUTPUT_DIR/Bdef.h" 2>/dev/null; then
        echo "  OK: QREAL typedef found"
    else
        echo "  MISSING: QREAL typedef"
        ((ERRORS++))
    fi

    if grep -q 'QCOMPLEX' "$OUTPUT_DIR/Bdef.h" 2>/dev/null; then
        echo "  OK: QCOMPLEX typedef found"
    else
        echo "  MISSING: QCOMPLEX typedef"
        ((ERRORS++))
    fi

    if grep -q 'MPI_QREAL' "$OUTPUT_DIR/Bdef.h" 2>/dev/null; then
        echo "  OK: MPI_QREAL declaration found"
    else
        echo "  MISSING: MPI_QREAL declaration"
        ((ERRORS++))
    fi

    echo ""
    echo "## Summary"
    echo "  Errors: $ERRORS"

} > "$VERIFY_REPORT" 2>&1

cat "$VERIFY_REPORT"
echo ""
echo "Report: $VERIFY_REPORT"
