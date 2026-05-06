#!/bin/bash
# Rename Fortran-mangled symbols (lowercase + trailing underscore) in
# a static archive from `<name>_` to `<name>_quad_`. Used by
# tests/blas/refblas and tests/lapack/reflapack to make refblas_quad /
# reflapack_quad coexist in the same link with the standard archive
# (which provides the un-renamed `<name>_` symbols).
#
# Renaming is restricted to a whitelist (one routine name per line,
# lowercase, no trailing underscore — written by refquad_alias.py
# alongside each rewritten ref_quad_*.f90). Routines outside the
# whitelist (typically precision-independent helpers — lsame_,
# xerbla_, ilaenv_, ieeeck_, etc.) keep their original names so the
# migrated test cycles, which don't include the standard archive in
# the link group, can still resolve those calls.
#
# objcopy --redefine-syms renames both definitions and references, so
# the archive stays self-consistent for the renamed routines'
# cross-calls.
#
# Usage: refquad_rename_archive.sh <archive.a> <whitelist.txt> [<more.txt>...]
#
# A stamp file alongside the archive marks completion so re-runs are
# no-ops (objcopy is non-idempotent without the guard — running twice
# yields `<name>_quad_quad_`).

set -e
ar=$1
shift
if [ -z "$ar" ] || [ ! -f "$ar" ]; then
    echo "usage: $0 <archive.a> <whitelist.txt> [<more.txt>...]" >&2
    exit 2
fi
if [ "$#" -lt 1 ]; then
    echo "$0: at least one whitelist file required" >&2
    exit 2
fi

stamp="${ar}.quad-renamed"
if [ -f "$stamp" ] && [ "$stamp" -nt "$ar" ]; then
    # Stamp newer than archive — already renamed since the last build.
    exit 0
fi

# Concatenate whitelists, deduplicate, build objcopy redefine list.
# Both definitions AND undefined references are renamed uniformly so
# cross-archive call sites stay consistent (e.g. reflapack_quad's
# `U dnrm2_` reference is renamed to `U dnrm2_quad_`, then resolved
# from refblas_quad which now exports `T dnrm2_quad_`).
whitelist=$(mktemp)
trap "rm -f $whitelist" EXIT
cat "$@" | tr -d '\r' | sed '/^\s*$/d' | sort -u \
    | awk '{print $1 "_ " $1 "_quad_"}' > "$whitelist"

if [ ! -s "$whitelist" ]; then
    echo "$0: empty whitelist (input: $*)" >&2
    touch "$stamp"
    exit 0
fi

n=$(wc -l < "$whitelist")
objcopy --redefine-syms="$whitelist" "$ar"
touch "$stamp"
echo "$0: renamed up to $n symbols in $(basename "$ar") (-> *_quad_)"
