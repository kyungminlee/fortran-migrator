#!/usr/bin/env bash
# 01-extract-symbols.sh — Extract and classify BLAS symbols from source.
#
# Scans BLAS source files for SUBROUTINE/FUNCTION definitions,
# classifies them by precision prefix, and writes a symbol report.
#
# This is a lightweight shell-based scanner for symbol extraction only.
# The fortran-migrator C++ tool has its own SymbolDatabase that does the
# same thing (and can also extract from compiled .a/.so libraries via nm).
# This script is primarily for diagnostics and report generation.
#
# For the full migration (type declarations, call sites, literals,
# intrinsics), the C++ source scanner is required — those patterns are
# too complex for reliable shell-based extraction.
#
# Usage: ./recipes/blas/01-extract-symbols.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

SYMBOL_REPORT="${REPORT_DIR}/symbols.txt"
RENAME_MAP="${REPORT_DIR}/rename-map.txt"

echo "=== Extracting BLAS symbols ==="

# --------------------------------------------------------------------------
# Extract SUBROUTINE and FUNCTION definitions from source
# --------------------------------------------------------------------------
extract_definitions() {
    local src_dir="$1"
    for f in "$src_dir"/*.f "$src_dir"/*.f90; do
        [[ -f "$f" ]] || continue
        local bname
        bname="$(basename "$f")"
        # Extract routine names from SUBROUTINE/FUNCTION definitions.
        # Skip comment lines (C/c/* in col 1 for .f, ! for .f90).
        # Use grep to find definition lines, then sed to extract
        # just the routine name (the identifier after SUBROUTINE or FUNCTION).
        grep -iE '(SUBROUTINE|FUNCTION)\s+[A-Za-z]' "$f" 2>/dev/null \
            | grep -vi '^[Cc*!]' \
            | grep -vi '^\s*!' \
            | sed -E 's/.*(SUBROUTINE|FUNCTION|subroutine|function|Subroutine|Function)[[:space:]]+([A-Za-z][A-Za-z0-9]*).*/\2/' \
            | tr 'a-z' 'A-Z' \
            | while read -r name; do
                echo "${name}	${bname}"
            done
    done
}

# --------------------------------------------------------------------------
# Classify precision prefix
# --------------------------------------------------------------------------
classify_symbol() {
    local name="$1"
    local first="${name:0:1}"
    local second="${name:1:1}"

    # ScaLAPACK: P + precision prefix
    if [[ "$first" == "P" ]] && [[ "$second" =~ [SDCZ] ]]; then
        echo "P${second}"
        return
    fi

    case "$first" in
        S) echo "S" ;;  # Single Real
        D) echo "D" ;;  # Double Real
        C) echo "C" ;;  # Single Complex
        Z) echo "Z" ;;  # Double Complex
        *) echo "-" ;;  # No precision prefix
    esac
}

# --------------------------------------------------------------------------
# Compute target prefix
# --------------------------------------------------------------------------
target_prefix() {
    local prefix="$1"
    local kind="$2"

    case "$prefix" in
        S|D)
            if [[ "$kind" == "10" ]]; then echo "E"; else echo "Q"; fi
            ;;
        C|Z)
            if [[ "$kind" == "10" ]]; then echo "Y"; else echo "X"; fi
            ;;
        PS|PD)
            if [[ "$kind" == "10" ]]; then echo "PE"; else echo "PQ"; fi
            ;;
        PC|PZ)
            if [[ "$kind" == "10" ]]; then echo "PY"; else echo "PX"; fi
            ;;
        *)
            echo "-"
            ;;
    esac
}

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
echo "Scanning: $BLAS_SRC"

{
    echo "# BLAS Symbol Database"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "# Source: $BLAS_SRC"
    echo "# Target KIND: $TARGET_KIND"
    echo "#"
    echo "# Format: NAME<tab>PREFIX<tab>BASE<tab>FILE"
    echo "#"

    extract_definitions "$BLAS_SRC" | sort -u | while IFS=$'\t' read -r name file; do
        prefix=$(classify_symbol "$name")
        if [[ "$prefix" == "-" ]]; then
            base="$name"
        elif [[ ${#prefix} -eq 2 ]]; then
            # ScaLAPACK P-prefix: base = P + rest after precision char
            base="P${name:2}"
        else
            base="${name:1}"
        fi
        echo -e "${name}\t${prefix}\t${base}\t${file}"
    done
} > "$SYMBOL_REPORT"

# --------------------------------------------------------------------------
# Generate rename map
# --------------------------------------------------------------------------
{
    echo "# BLAS Rename Map (KIND=$TARGET_KIND)"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "#"
    echo "# Format: ORIGINAL<tab>TARGET<tab>PREFIX<tab>FILE"
    echo "#"

    grep -v '^#' "$SYMBOL_REPORT" | while IFS=$'\t' read -r name prefix base file; do
        tp=$(target_prefix "$prefix" "$TARGET_KIND")
        if [[ "$tp" == "-" ]]; then
            target="$name"
        elif [[ ${#tp} -eq 2 ]]; then
            target="${tp}${base:1}"
        else
            target="${tp}${base}"
        fi
        echo -e "${name}\t${target}\t${prefix}\t${file}"
    done
} > "$RENAME_MAP"

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
total=$(grep -vc '^#' "$SYMBOL_REPORT" || true)
single_r=$(grep -c $'\tS\t' "$SYMBOL_REPORT" || true)
double_r=$(grep -c $'\tD\t' "$SYMBOL_REPORT" || true)
single_c=$(grep -c $'\tC\t' "$SYMBOL_REPORT" || true)
double_c=$(grep -c $'\tZ\t' "$SYMBOL_REPORT" || true)
none=$(grep -c $'\t-\t' "$SYMBOL_REPORT" || true)

echo ""
echo "Symbols extracted: $total"
echo "  Single Real (S):    $single_r"
echo "  Double Real (D):    $double_r"
echo "  Single Complex (C): $single_c"
echo "  Double Complex (Z): $double_c"
echo "  No prefix:          $none"
echo ""
echo "Symbol report: $SYMBOL_REPORT"
echo "Rename map:    $RENAME_MAP"
