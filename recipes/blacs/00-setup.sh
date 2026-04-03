#!/usr/bin/env bash
# 00-setup.sh — Set up directories and variables for BLACS migration.
#
# Usage: source recipes/blacs/00-setup.sh [--kind 10|16]

set -euo pipefail

RECIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${RECIPE_DIR}/../.." && pwd)"

TARGET_KIND="${TARGET_KIND:-16}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --kind) TARGET_KIND="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; return 1 2>/dev/null || exit 1 ;;
    esac
done

if [[ "$TARGET_KIND" != "10" && "$TARGET_KIND" != "16" ]]; then
    echo "Error: --kind must be 10 or 16" >&2
    return 1 2>/dev/null || exit 1
fi

# Prefix mapping
if [[ "$TARGET_KIND" == "10" ]]; then
    REAL_PREFIX="e"
    COMPLEX_PREFIX="y"
    REAL_PREFIX_UPPER="E"
    COMPLEX_PREFIX_UPPER="Y"
else
    REAL_PREFIX="q"
    COMPLEX_PREFIX="x"
    REAL_PREFIX_UPPER="Q"
    COMPLEX_PREFIX_UPPER="X"
fi

BLACS_SRC="${PROJECT_ROOT}/external/scalapack-2.2.3/BLACS/SRC"
WORK_DIR="${RECIPE_DIR}/work"
OUTPUT_DIR="${WORK_DIR}/output-kind${TARGET_KIND}"
REPORT_DIR="${WORK_DIR}/reports-kind${TARGET_KIND}"
BUILD_DIR="${WORK_DIR}/build-kind${TARGET_KIND}"

if [[ ! -d "$BLACS_SRC" ]]; then
    echo "Error: BLACS source not found: $BLACS_SRC" >&2
    return 1 2>/dev/null || exit 1
fi

mkdir -p "$OUTPUT_DIR" "$REPORT_DIR" "$BUILD_DIR"

echo "BLACS migration setup:"
echo "  TARGET_KIND      = $TARGET_KIND"
echo "  REAL_PREFIX       = $REAL_PREFIX (${REAL_PREFIX_UPPER})"
echo "  COMPLEX_PREFIX    = $COMPLEX_PREFIX (${COMPLEX_PREFIX_UPPER})"
echo "  BLACS_SRC         = $BLACS_SRC"
echo "  OUTPUT_DIR        = $OUTPUT_DIR"

export TARGET_KIND REAL_PREFIX COMPLEX_PREFIX REAL_PREFIX_UPPER COMPLEX_PREFIX_UPPER
export BLACS_SRC OUTPUT_DIR REPORT_DIR BUILD_DIR PROJECT_ROOT RECIPE_DIR WORK_DIR
