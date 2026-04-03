#!/usr/bin/env bash
# 02-migrate.sh — Migrate BLAS source files to extended precision.
#
# Uses the Python migration engine (migrate_blas.py) to convert all
# BLAS source files to the target KIND precision.
#
# Usage: ./recipes/blas/02-migrate.sh [--kind 10|16] [--dry-run]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DRY_RUN=""
KIND_ARGS=()

for arg in "$@"; do
    if [[ "$arg" == "--dry-run" ]]; then
        DRY_RUN="--dry-run"
    elif [[ "$arg" == "--kind" ]]; then
        KIND_ARGS+=("$arg")
    elif [[ "${KIND_ARGS[-1]:-}" == "--kind" ]]; then
        KIND_ARGS+=("$arg")
    fi
done

source "${SCRIPT_DIR}/00-setup.sh" "${KIND_ARGS[@]}"

echo "=== Migrating BLAS to KIND=$TARGET_KIND ==="

# Clean output directory
rm -rf "${OUTPUT_DIR:?}"/*

echo "Source:  $BLAS_SRC"
echo "Output:  $OUTPUT_DIR"
echo ""

# Run the Python migration engine
python3 "${SCRIPT_DIR}/migrate_blas.py" \
    --kind "$TARGET_KIND" \
    $DRY_RUN \
    "$BLAS_SRC" \
    "$OUTPUT_DIR" \
    2>&1 | tee "${REPORT_DIR}/migrate.log"

# --------------------------------------------------------------------------
# Post-migration summary
# --------------------------------------------------------------------------
if [[ -z "$DRY_RUN" ]]; then
    echo ""
    echo "=== Migration Output ==="
    input_count=$(find "$BLAS_SRC" -maxdepth 1 \( -name '*.f' -o -name '*.f90' \) | wc -l | tr -d ' ')
    output_count=$(find "$OUTPUT_DIR" -maxdepth 1 \( -name '*.f' -o -name '*.f90' \) 2>/dev/null | wc -l | tr -d ' ')
    echo "Input files:  $input_count"
    echo "Output files: $output_count"

    # List renamed files
    if [[ "$output_count" -gt 0 ]]; then
        echo ""
        echo "Output files:"
        ls "$OUTPUT_DIR"/*.f "$OUTPUT_DIR"/*.f90 2>/dev/null | xargs -I{} basename {} | sort | column
    fi
fi
