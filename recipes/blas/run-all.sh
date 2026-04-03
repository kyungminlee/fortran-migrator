#!/usr/bin/env bash
# run-all.sh — Run the full BLAS migration pipeline.
#
# Executes all steps in order:
#   01 - Extract symbols
#   02 - Migrate source files
#   03 - Convergence check
#   04 - Verify migrated output
#   05 - Compile
#   06 - Build library
#
# Usage: ./recipes/blas/run-all.sh [--kind 10|16]
#
# Set SKIP_COMPILE=1 to skip compilation steps (04, 05, 06).
# Set STOP_ON_ERROR=0 to continue past failures.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STOP_ON_ERROR="${STOP_ON_ERROR:-1}"
SKIP_COMPILE="${SKIP_COMPILE:-0}"

run_step() {
    local script="$1"
    shift
    local name
    name=$(basename "$script" .sh)
    echo ""
    echo "================================================================"
    echo "  Step: $name"
    echo "================================================================"
    echo ""

    if bash "$script" "$@"; then
        echo ""
        echo ">>> $name: OK"
    else
        echo ""
        echo ">>> $name: FAILED"
        if [[ "$STOP_ON_ERROR" == "1" ]]; then
            echo "Stopping. Set STOP_ON_ERROR=0 to continue past failures."
            exit 1
        fi
    fi
}

echo "========================================"
echo " BLAS Migration Pipeline"
echo "========================================"
echo " Arguments: $*"

run_step "${SCRIPT_DIR}/01-extract-symbols.sh" "$@"
run_step "${SCRIPT_DIR}/02-migrate.sh" "$@"
run_step "${SCRIPT_DIR}/03-convergence-check.sh" "$@"
run_step "${SCRIPT_DIR}/04-verify.sh" "$@"

if [[ "$SKIP_COMPILE" != "1" ]]; then
    run_step "${SCRIPT_DIR}/05-compile.sh" "$@"
    run_step "${SCRIPT_DIR}/06-build-library.sh" "$@"
else
    echo ""
    echo "Skipping compilation steps (SKIP_COMPILE=1)"
fi

echo ""
echo "========================================"
echo " Pipeline complete."
echo "========================================"
