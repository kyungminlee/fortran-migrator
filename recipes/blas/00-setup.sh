#!/usr/bin/env bash
# 00-setup.sh — Set up directories and variables for BLAS migration.
#
# Usage: source recipes/blas/00-setup.sh [--kind 10|16]
#
# This script is meant to be sourced (not executed) so that the
# environment variables it defines are available to subsequent scripts.

set -euo pipefail

# --------------------------------------------------------------------------
# Project root (parent of recipes/)
# --------------------------------------------------------------------------
RECIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${RECIPE_DIR}/../.." && pwd)"

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
TARGET_KIND="${TARGET_KIND:-16}"

# Parse optional arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --kind)
            TARGET_KIND="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            return 1 2>/dev/null || exit 1
            ;;
    esac
done

if [[ "$TARGET_KIND" != "10" && "$TARGET_KIND" != "16" ]]; then
    echo "Error: --kind must be 10 or 16 (got: $TARGET_KIND)" >&2
    return 1 2>/dev/null || exit 1
fi

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
BLAS_SRC="${PROJECT_ROOT}/external/lapack-3.12.1/BLAS/SRC"
MIGRATOR="${PROJECT_ROOT}/build/fortran-migrator"

# Output directories
WORK_DIR="${RECIPE_DIR}/work"
OUTPUT_DIR="${WORK_DIR}/output-kind${TARGET_KIND}"
CONVERGENCE_DIR="${WORK_DIR}/convergence-kind${TARGET_KIND}"
REPORT_DIR="${WORK_DIR}/reports-kind${TARGET_KIND}"
BUILD_DIR="${WORK_DIR}/build-kind${TARGET_KIND}"

# --------------------------------------------------------------------------
# Validate
# --------------------------------------------------------------------------
if [[ ! -d "$BLAS_SRC" ]]; then
    echo "Error: BLAS source directory not found: $BLAS_SRC" >&2
    echo "  Run: cd ${PROJECT_ROOT}/external && tar xzf lapack-3.12.1.tar.gz" >&2
    return 1 2>/dev/null || exit 1
fi

if [[ ! -x "$MIGRATOR" ]]; then
    echo "Warning: fortran-migrator not built: $MIGRATOR" >&2
    echo "  Run: cd ${PROJECT_ROOT}/build && cmake --build ." >&2
fi

# --------------------------------------------------------------------------
# Create directories
# --------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR" "$CONVERGENCE_DIR" "$REPORT_DIR" "$BUILD_DIR"

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
echo "BLAS migration setup:"
echo "  TARGET_KIND    = $TARGET_KIND"
echo "  BLAS_SRC       = $BLAS_SRC"
echo "  OUTPUT_DIR     = $OUTPUT_DIR"
echo "  CONVERGENCE_DIR= $CONVERGENCE_DIR"
echo "  REPORT_DIR     = $REPORT_DIR"
echo "  BUILD_DIR      = $BUILD_DIR"
echo "  MIGRATOR       = $MIGRATOR"

# Export for use by other scripts
export TARGET_KIND BLAS_SRC MIGRATOR OUTPUT_DIR CONVERGENCE_DIR REPORT_DIR BUILD_DIR PROJECT_ROOT RECIPE_DIR WORK_DIR
