#!/usr/bin/env bash
# 01-extract-symbols.sh — Extract and classify BLACS C symbols.
#
# Scans BLACS C source files for function definitions,
# classifies them by precision prefix and role.
#
# Usage: ./recipes/blacs/01-extract-symbols.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

SYMBOL_REPORT="${REPORT_DIR}/symbols.txt"

echo "=== Extracting BLACS symbols ==="

{
    echo "# BLACS Symbol Database"
    echo "# Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "#"
    echo "# Category: user=user-facing, vv=vector-vector, mpi=MPI op, copy=matrix copy, infra=infrastructure"
    echo "# Format: FILE<tab>PREFIX<tab>CATEGORY<tab>FUNCTION_NAME"
    echo "#"

    cd "$BLACS_SRC"

    # User-facing routines: [sdczi]ge{sd,rv,bs,br}2d_, [sdczi]tr{sd,rv,bs,br}2d_, [sdczi]g{sum,amx,amn}2d_
    for f in [sdczi]ge*2d_.c [sdczi]tr*2d_.c [sdczi]g*2d_.c; do
        [[ -f "$f" ]] || continue
        prefix="${f:0:1}"
        case "$prefix" in
            s|d|c|z|i) ;;
            *) continue ;;
        esac
        func=$(basename "$f" .c)
        echo -e "${f}\t${prefix}\tuser\t${func}"
    done

    # Vector-vector operations: BI_[sdczi]vv{sum,amn,amn2,amx,amx2}.c
    for f in BI_[sdczi]vv*.c; do
        [[ -f "$f" ]] || continue
        prefix="${f:3:1}"
        func=$(grep -oE 'void BI_[a-z]+[a-z0-9]*' "$f" 2>/dev/null | head -1 | sed 's/void //' || echo "")
        echo -e "${f}\t${prefix}\tvv\t${func}"
    done

    # MPI operation wrappers: BI_[sdczi]MPI_*.c
    for f in BI_[sdczi]MPI_*.c; do
        [[ -f "$f" ]] || continue
        prefix="${f:3:1}"
        func=$(grep -oE 'void BI_[a-zA-Z_]+[a-z0-9]*' "$f" 2>/dev/null | head -1 | sed 's/void //' || echo "")
        echo -e "${f}\t${prefix}\tmpi\t${func}"
    done

    # Matrix copy: BI_[sd]mvcopy.c, BI_[sd]vmcopy.c
    for f in BI_[sd]mvcopy.c BI_[sd]vmcopy.c; do
        [[ -f "$f" ]] || continue
        prefix="${f:3:1}"
        func=$(grep -oE 'void BI_[a-z]+' "$f" 2>/dev/null | head -1 | sed 's/void //' || echo "")
        echo -e "${f}\t${prefix}\tcopy\t${func}"
    done

    # Infrastructure (type-independent)
    for f in BI_[A-Z]*.c BI_i*.c blacs_*.c sys2blacs*.c; do
        [[ -f "$f" ]] || continue
        # Skip already-classified files
        [[ "$f" == BI_[sdcz]vv*.c ]] && continue
        [[ "$f" == BI_[sdcz]MPI_*.c ]] && continue
        [[ "$f" == BI_[sd]mvcopy.c ]] && continue
        [[ "$f" == BI_[sd]vmcopy.c ]] && continue
        echo -e "${f}\t-\tinfra\t-"
    done

} | sort > "$SYMBOL_REPORT"

# Summary
total=$(grep -vc '^#' "$SYMBOL_REPORT" || true)
user=$(grep -c $'\tuser\t' "$SYMBOL_REPORT" || true)
vv=$(grep -c $'\tvv\t' "$SYMBOL_REPORT" || true)
mpi_op=$(grep -c $'\tmpi\t' "$SYMBOL_REPORT" || true)
copy=$(grep -c $'\tcopy\t' "$SYMBOL_REPORT" || true)
infra=$(grep -c $'\tinfra\t' "$SYMBOL_REPORT" || true)

echo ""
echo "Files classified: $total"
echo "  User-facing (s/d/c/z/i): $user"
echo "  Vector-vector ops:       $vv"
echo "  MPI op wrappers:         $mpi_op"
echo "  Matrix copy:             $copy"
echo "  Infrastructure:          $infra"
echo ""
echo "Files to clone (d→${REAL_PREFIX}, z→${COMPLEX_PREFIX}):"
d_count=$(grep $'\td\t' "$SYMBOL_REPORT" | grep -cv 'infra' || true)
z_count=$(grep $'\tz\t' "$SYMBOL_REPORT" | grep -cv 'infra' || true)
echo "  d-variants: $d_count → ${REAL_PREFIX}-variants"
echo "  z-variants: $z_count → ${COMPLEX_PREFIX}-variants"
echo ""
echo "Report: $SYMBOL_REPORT"
