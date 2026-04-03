#!/usr/bin/env bash
# 02-migrate.sh — Migrate BLACS C source files to extended precision.
#
# Uses the general-purpose pyengine pipeline with recipes/blacs.yaml.
#
# Usage: ./recipes/blacs/02-migrate.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Migrating BLACS to KIND=$TARGET_KIND ==="

# Clean output
rm -rf "${OUTPUT_DIR:?}"/*

echo "Source:  $BLACS_SRC"
echo "Output:  $OUTPUT_DIR"
echo ""

# Run the general-purpose migration engine
uv run python -m pyengine \
    "${PROJECT_ROOT}/recipes/blacs.yaml" \
    "$OUTPUT_DIR" \
    --kind "$TARGET_KIND" \
    --project-root "$PROJECT_ROOT" \
    2>&1 | tee "${REPORT_DIR}/migrate.log"

echo ""
orig_count=$(find "$BLACS_SRC" -maxdepth 1 -name '*.c' | wc -l | tr -d ' ')
new_count=$(find "$OUTPUT_DIR" -maxdepth 1 -name '*.c' 2>/dev/null | wc -l | tr -d ' ')
added=$((new_count - orig_count))
echo "Original files: $orig_count"
echo "Output files:   $new_count (+$added new)"
