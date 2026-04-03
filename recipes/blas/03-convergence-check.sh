#!/usr/bin/env bash
# 03-convergence-check.sh — Verify S→target and D→target produce identical output.
#
# Migrates from both S and D (and C and Z) variants to the same target KIND,
# then diffs each convergence pair. Differences reveal:
#   - Comment-only: harmless (different wording in doc headers)
#   - Algorithm: intentional differences between S/D implementations
#   - Missed conversions: bugs in the migrator
#
# Usage: ./recipes/blas/03-convergence-check.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Convergence Check (KIND=$TARGET_KIND) ==="

# --------------------------------------------------------------------------
# Determine target prefixes
# --------------------------------------------------------------------------
if [[ "$TARGET_KIND" == "10" ]]; then
    REAL_PREFIX="e"
    COMPLEX_PREFIX="y"
else
    REAL_PREFIX="q"
    COMPLEX_PREFIX="x"
fi

# --------------------------------------------------------------------------
# Define convergence pairs
# --------------------------------------------------------------------------
# Format: base_name:s_file:d_file:target_file
# Real pairs (S↔D → Q/E)
REAL_PAIRS=(
    asum axpy copy dot gbmv gemm gemmtr gemv ger
    sbmv scal spmv spr spr2 swap symm symv syr syr2
    syr2k syrk tbmv tbsv tpmv tpsv trmm trmv trsm trsv
    rot rotm rotmg
)

# Complex pairs (C↔Z → X/Y)
COMPLEX_PAIRS=(
    axpy copy dotc dotu gbmv gemm gemmtr gemv gerc geru
    hbmv hemm hemv her her2 her2k herk hpmv hpr hpr2
    scal swap symm syr2k syrk tbmv tbsv tpmv tpsv
    trmm trmv trsm trsv
)

# F90 real pairs
F90_REAL_PAIRS=(nrm2)

# F90 complex pairs
F90_COMPLEX_PAIRS=(rotg)

# Cross-type pairs
# Format: s_file:d_file:target_file
CROSS_PAIRS=(
    "scabs1:dcabs1:${REAL_PREFIX}cabs1"
    "scasum:dzasum:${REAL_PREFIX}${COMPLEX_PREFIX}asum"
    "isamax:idamax:i${REAL_PREFIX}amax"
    "icamax:izamax:i${COMPLEX_PREFIX}amax"
)

# F90 cross-type pairs
F90_CROSS_PAIRS=(
    "scnrm2:dznrm2:${REAL_PREFIX}${COMPLEX_PREFIX}nrm2"
)

# --------------------------------------------------------------------------
# Run convergence diffs
# --------------------------------------------------------------------------
CONV_REPORT="${REPORT_DIR}/convergence.txt"
PASS=0
FAIL=0
COMMENT_ONLY=0
MISSING=0

{
    echo "# Convergence Check Report"
    echo "# TARGET_KIND=$TARGET_KIND"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo ""

    # --- Real pairs (.f) ---
    for base in "${REAL_PAIRS[@]}"; do
        s_file="${OUTPUT_DIR}/${REAL_PREFIX}${base}.f"
        d_file="${OUTPUT_DIR}/${REAL_PREFIX}${base}.f"  # Both S and D map to same target

        # The migrated S-variant and D-variant should both exist as the same target name.
        # To compare, we need separate migration output directories.
        # For now, check if the output file exists.
        target="${REAL_PREFIX}${base}.f"

        if [[ ! -f "${OUTPUT_DIR}/${target}" ]]; then
            echo "MISSING: ${target} (from s${base}.f / d${base}.f)"
            ((MISSING++))
            continue
        fi

        # If we had separate S-origin and D-origin output, diff them here.
        # For now, mark as needing --convergence-check mode in the migrator.
        echo "PENDING: ${target} (requires --convergence-check mode)"
    done

    # --- Real pairs (.f90) ---
    for base in "${F90_REAL_PAIRS[@]}"; do
        target="${REAL_PREFIX}${base}.f90"
        if [[ ! -f "${OUTPUT_DIR}/${target}" ]]; then
            echo "MISSING: ${target}"
            ((MISSING++))
        else
            echo "PENDING: ${target} (requires --convergence-check mode)"
        fi
    done

    # --- Complex pairs (.f) ---
    for base in "${COMPLEX_PAIRS[@]}"; do
        target="${COMPLEX_PREFIX}${base}.f"
        if [[ ! -f "${OUTPUT_DIR}/${target}" ]]; then
            echo "MISSING: ${target} (from c${base}.f / z${base}.f)"
            ((MISSING++))
        else
            echo "PENDING: ${target} (requires --convergence-check mode)"
        fi
    done

    # --- Complex pairs (.f90) ---
    for base in "${F90_COMPLEX_PAIRS[@]}"; do
        target="${COMPLEX_PREFIX}${base}.f90"
        if [[ ! -f "${OUTPUT_DIR}/${target}" ]]; then
            echo "MISSING: ${target}"
            ((MISSING++))
        else
            echo "PENDING: ${target} (requires --convergence-check mode)"
        fi
    done

    echo ""
    echo "# Summary"
    echo "# PASS:         $PASS"
    echo "# FAIL:         $FAIL"
    echo "# COMMENT_ONLY: $COMMENT_ONLY"
    echo "# MISSING:      $MISSING"

} > "$CONV_REPORT"

echo ""
cat "$CONV_REPORT"
echo ""
echo "Report: $CONV_REPORT"
echo ""
echo "Note: Full convergence checking requires the migrator's"
echo "  --convergence-check mode, which migrates from both S and D"
echo "  (or C and Z) independently and diffs the results."
