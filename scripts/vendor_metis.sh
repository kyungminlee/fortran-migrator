#!/bin/bash
#
# vendor_metis.sh — regenerate extern/metis-5.2.1 from upstream.
#
# METIS is vendored rather than found on the system so MUMPS's ICNTL(7)=5
# nested-dissection ordering is always available, and so the copy we link
# can be *privately namespaced*: every linker-visible symbol carries a
# MUMPS-specific prefix, and the archive therefore cannot clash (ODR /
# duplicate-symbol) with a system METIS that the final application might
# also pull in.
#
# The vendored tree is checked into git already renamed. This script is
# how that tree is produced, so the rename stays mechanical, auditable,
# and repeatable when upstream moves. It is a maintainer tool — the build
# never runs it.
#
# Usage:
#   scripts/vendor_metis.sh [<dest>]
#
#   <dest>  output directory (default: extern/metis-5.2.1). Removed and
#           recreated. Nothing else in the repo is touched — rewiring
#           cmake/EplinalgMumps.cmake and codegen/migrator/staging.py to
#           a new directory name is a manual follow-up.
#
# Environment:
#   CC       compiler used to enumerate GKlib's exports and to verify the
#            result (default: cc). Required — see step 5; the rename it
#            drives is part of producing the tree, not just checking it.
#   WORKDIR  scratch directory for the upstream clones (default: a mktemp
#            dir, removed on exit). Point it at a persistent path to
#            avoid re-cloning while iterating.
#
# Upstream layout, from 5.2.1 on, is two repositories: METIS no longer
# bundles GKlib. Both are pinned below — GKlib has no release tags at
# all, so it is pinned by commit sha.
#
# The rename, in full:
#
#   1. Public API.  The 13 entry points declared METIS_API in
#      include/metis.h, in all the casings the Fortran wrappers emit:
#        METIS_NodeND   -> METIS_MUMPS_NodeND
#        METIS_NODEND   -> METIS_MUMPS_NODEND
#        metis_nodend{,_,__} -> metis_mumps_nodend{,_,__}
#      The stem list is read out of metis.h, not hardcoded, so a new
#      upstream entry point is picked up automatically.
#
#   2. libmetis internals.  Upstream already funnels them through
#      libmetis/rename.h and libmetis/gklib_rename.h with a libmetis__
#      prefix; that prefix is rewritten to libmetis_MUMPS_.
#
#   3. Internals upstream forgot to list in rename.h (LIBMETIS_STRAYS).
#      Appended as extra #defines rather than edited in place, matching
#      how upstream handles the rest.
#
#   4. Every GKlib global — the ~640 that carry GKlib's own gk_ prefix
#      plus the couple dozen that carry none (GK_STRAYS) — gets MUMPS_
#      *prepended*: gk_malloc becomes MUMPS_gk_malloc, errexit becomes
#      MUMPS_errexit. gk_ is GKlib's prefix, not ours: any other GKlib in
#      the link — a system one, ParMETIS's, METIS's own programs —
#      exports the same names. Prepending rather than substituting the
#      gk_ keeps the mapping injective; GKlib ships both errexit and
#      gk_errexit, two live functions with different signatures, and
#      substitution collapses them onto one name. The gk_ list is read
#      out of a throwaway build rather than hardcoded, so it cannot rot,
#      and only actual defined globals are renamed — gk_ types, macros,
#      and struct members are untouched.
#
#   5. metis.h's IDXTYPEWIDTH / REALTYPEWIDTH self-defaults, which 5.2.1
#      comments out, are restored to 32 — the width this project builds
#      and the 5.1.0 behavior. Without them every translation unit that
#      includes the *installed* metis.h would have to supply -D itself.
#
#   6. GKlib's io.c defines _GNU_SOURCE around its <stdio.h> include,
#      purely to get getline(). Since glibc 2.38 that same macro also
#      switches on the C23 scanf redirects, so io.c alone in the whole
#      tree emits a call to __isoc23_sscanf — a symbol that simply does
#      not exist in older glibc, which makes the archive unlinkable
#      against them (this project targets a glibc 2.28 floor). getline
#      is POSIX-2008, so _POSIX_C_SOURCE=200809L declares it without
#      dragging in the GNU extensions, and sscanf stays __isoc99_sscanf
#      (glibc 2.7+).
#
# Everything else is upstream verbatim; `diff -r` against a pristine
# checkout with the rename reversed comes back empty.
#
# Verification: the result is compiled with this project's flags and
# every one of the ~980 defined globals in the resulting objects is
# checked to carry one of the four private prefixes — METIS_MUMPS_,
# metis_mumps_, libmetis_MUMPS_, MUMPS_. A brand-new unprefixed name
# aborts the run and names itself — add it to LIBMETIS_STRAYS or
# GK_STRAYS and re-run.
#
set -euo pipefail

METIS_REPO=https://github.com/KarypisLab/METIS.git
METIS_REF=v5.2.1
GKLIB_REPO=https://github.com/KarypisLab/GKlib.git
GKLIB_SHA=3b7d61b9f885063c89901f3901fb4426f9cfb58f

# Internals defined in libmetis/ that upstream's rename.h does not list.
LIBMETIS_STRAYS=(
    BalanceAndRefineLP
    BlockKWayPartitioning
    CoarsenGraphNlevels
    ComputeBFSOrdering
    Greedy_KWayEdgeCutOptimize
    Greedy_KWayEdgeStats
    GrowBisectionNode2
    GrowMultisection
)

# GKlib exports that lack GKlib's own gk_ prefix.
GK_STRAYS=(
    ComputeAccuracy ComputeMean ComputeMedianRFP ComputeROCn ComputeStdDev
    decodeblock encodeblock errexit getpathname
    GKDecodeBase64 GKEncodeBase64 gkfooo
    HTable_Create HTable_Delete HTable_Destroy HTable_GetNext
    HTable_HFunction HTable_Insert HTable_Reset HTable_Resize
    HTable_Search HTable_SearchAndDelete
    itemsets_find_frequent_itemsets itemsets_project_matrix
    PrintBackTrace
)

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
dest=${1:-$root/extern/metis-5.2.1}
case $dest in /*) ;; *) dest=$root/$dest ;; esac

if [[ -n ${WORKDIR:-} ]]; then
    work=$WORKDIR
    mkdir -p "$work"
else
    work=$(mktemp -d)
    trap 'rm -rf "$work"' EXIT
fi

say() { printf '  %s\n' "$*"; }

# ── 1. fetch upstream ───────────────────────────────────────────────
if [[ ! -d $work/METIS/.git ]]; then
    say "cloning METIS $METIS_REF"
    git -c advice.detachedHead=false clone --quiet --depth 1 --branch "$METIS_REF" "$METIS_REPO" "$work/METIS"
fi
if [[ ! -d $work/GKlib/.git ]]; then
    say "cloning GKlib $GKLIB_SHA"
    git clone --quiet "$GKLIB_REPO" "$work/GKlib"
    git -C "$work/GKlib" checkout --quiet "$GKLIB_SHA"
fi
have=$(git -C "$work/GKlib" rev-parse HEAD)
[[ $have == "$GKLIB_SHA" ]] || { echo "GKlib at $have, expected $GKLIB_SHA" >&2; exit 1; }

# ── 2. assemble the pruned tree ─────────────────────────────────────
# Kept: the sources this project compiles, the public header, and the
# licensing / provenance files. Dropped: build systems (we supply our
# own), the standalone programs, tests, sample graphs, and the manual.
say "assembling $dest"
rm -rf "$dest"
mkdir -p "$dest/include" "$dest/libmetis" "$dest/GKlib/include" "$dest/GKlib/src"

cp "$work/METIS/include/metis.h"          "$dest/include/"
cp "$work/METIS"/libmetis/*.c             "$dest/libmetis/"
cp "$work/METIS"/libmetis/*.h             "$dest/libmetis/"
cp "$work/METIS"/{LICENSE,Changelog,README.md} "$dest/"

cp "$work/GKlib"/include/*.h              "$dest/GKlib/include/"
cp "$work/GKlib"/src/*.c                  "$dest/GKlib/src/"
cp "$work/GKlib"/{LICENSE.txt,LICENSES.md,README.md} "$dest/GKlib/"
cp -r "$work/GKlib/LICENSES"              "$dest/GKlib/"

# ── 3. local patches ────────────────────────────────────────────────
# 3a. restore the width self-defaults
sed -i -e 's|^//#define IDXTYPEWIDTH 32$|#define IDXTYPEWIDTH 32|' \
       -e 's|^//#define REALTYPEWIDTH 32$|#define REALTYPEWIDTH 32|' \
       "$dest/include/metis.h"
grep -q '^#define IDXTYPEWIDTH 32$'  "$dest/include/metis.h" || { echo "IDXTYPEWIDTH patch failed" >&2; exit 1; }
grep -q '^#define REALTYPEWIDTH 32$' "$dest/include/metis.h" || { echo "REALTYPEWIDTH patch failed" >&2; exit 1; }

# 3b. keep glibc's C23 scanf redirects out of io.c. Upstream reaches for
# _GNU_SOURCE only to declare getline(); since glibc 2.38 that macro also
# redirects sscanf to __isoc23_sscanf, which does not exist before 2.38
# and so pins the archive to a build host at least that new. getline is
# POSIX-2008 — ask for exactly that instead.
sed -i -e 's|^/\* Get getline to be defined\. \*/$|/* Get getline to be defined. _POSIX_C_SOURCE, not _GNU_SOURCE: the\n   latter also turns on glibc 2.38+ C23 scanf redirects (__isoc23_sscanf),\n   which do not exist in older glibc. See scripts/vendor_metis.sh. */|' \
       -e 's|^#define _GNU_SOURCE$|#define _POSIX_C_SOURCE 200809L|' \
       -e 's|^#undef _GNU_SOURCE$|#undef _POSIX_C_SOURCE|' \
       "$dest/GKlib/src/io.c"
grep -q '^#define _POSIX_C_SOURCE 200809L$' "$dest/GKlib/src/io.c" || { echo "io.c _GNU_SOURCE patch failed" >&2; exit 1; }
grep -qE '^#(define|undef) _GNU_SOURCE$' "$dest/GKlib/src/io.c" && { echo "io.c still defines _GNU_SOURCE" >&2; exit 1; }

# ── 4. rename ───────────────────────────────────────────────────────
# Public API stems, straight out of the header.
mapfile -t stems < <(grep -oE '^METIS_API\([a-z]+\) METIS_[A-Za-z]+' "$dest/include/metis.h" \
                     | sed 's/.*METIS_//' | sort -u)
[[ ${#stems[@]} -gt 0 ]] || { echo "no METIS_API entry points found in metis.h" >&2; exit 1; }
say "public API: ${#stems[@]} entry points"

sedprog=$work/rename.sed
: > "$sedprog"
for s in "${stems[@]}"; do
    up=${s^^}
    lo=${s,,}
    printf 's/\\bMETIS_%s\\b/METIS_MUMPS_%s/g\n'   "$s"  "$s"  >> "$sedprog"
    printf 's/\\bMETIS_%s\\b/METIS_MUMPS_%s/g\n'   "$up" "$up" >> "$sedprog"
    printf 's/\\bmetis_%s\\b/metis_mumps_%s/g\n'   "$lo" "$lo" >> "$sedprog"
    printf 's/\\bmetis_%s_\\b/metis_mumps_%s_/g\n' "$lo" "$lo" >> "$sedprog"
    printf 's/\\bmetis_%s__\\b/metis_mumps_%s__/g\n' "$lo" "$lo" >> "$sedprog"
done
# Internals upstream already routes through rename.h / gklib_rename.h.
printf 's/\\blibmetis__/libmetis_MUMPS_/g\n' >> "$sedprog"

find "$dest" \( -name '*.c' -o -name '*.h' \) -print0 | xargs -0 sed -i -f "$sedprog"

# Strays upstream never listed. Appended as extra #defines inside the
# existing include guards, so they behave exactly like the entries
# upstream did write.
{
    printf '\n/* Internals that escaped the list above (see scripts/vendor_metis.sh). */\n'
    for n in "${LIBMETIS_STRAYS[@]}"; do
        printf '#define %-32s libmetis_MUMPS_%s\n' "$n" "$n"
    done
} > "$work/libmetis_strays.h"
awk -v strays="$(cat "$work/libmetis_strays.h")" '
    /^#endif/ && !done { print strays; done = 1 }
    { print }
' "$dest/libmetis/rename.h" > "$work/rename.h.new"
mv "$work/rename.h.new" "$dest/libmetis/rename.h"

# gkregex.c defines its no-op stub *above* its own #include "GKlib.h",
# so the macro above cannot reach it; rename the definition in place.
sed -i 's/^void gkfooo()/void MUMPS_gkfooo()/' "$dest/GKlib/src/gkregex.c"
grep -q '^void MUMPS_gkfooo()' "$dest/GKlib/src/gkregex.c" || { echo "gkfooo patch failed" >&2; exit 1; }

# ── 5. namespace GKlib ──────────────────────────────────────────────
# One rule for every GKlib global: MUMPS_ is *prepended* to the name
# upstream chose. gk_malloc becomes MUMPS_gk_malloc; errexit, which never
# had GKlib's prefix, becomes MUMPS_errexit.
#
# Prepending rather than substituting the gk_ matters: GKlib ships both
# errexit(char *, ...) and gk_errexit(int signum, char *, ...) — two live
# functions with different signatures — and s/^gk_/MUMPS_/ collapses them
# onto one name, which is a compile error at best. Prepending is
# injective for free, keeps upstream's name legible inside ours, and
# reverses by stripping a fixed string.
#
# Steps 1-3 left the ~640 gk_-prefixed globals alone: gk_ is at least not
# *bare*. But gk_ is GKlib's prefix, not ours — a system GKlib, or any
# other library that vendors one (ParMETIS, METIS's own programs),
# exports the same names, and two gk_mallocs in one link is exactly the
# clash this script exists to prevent.
#
# Those ~640 are not hardcoded; a list that size would rot the first time
# upstream moved. It is read out of a throwaway build of the tree as it
# stands after step 4, which is why a compiler is required to produce the
# tree and not merely to check it. Only actual defined globals are
# renamed, so gk_ types, macros, and struct members are untouched — a
# blanket s/gk_/MUMPS_gk_/ would have churned all of those for nothing.
# GK_STRAYS, the globals with no prefix at all, stays an explicit list:
# the verify step rejects anything bare, and the fix is to name it there.
#
# The block lands before GKlib.h's own gk_*.h includes — not before
# <gk_proto.h>, where the strays alone used to go. The eight gk_ globals
# that are data rather than functions (gk_optind and friends) are
# declared in gk_getopt.h and gk_externs.h, and a #define placed after
# those headers renames the definitions while leaving the declarations
# behind. Every system header GKlib.h pulls in is already above that
# line, so nothing outside GKlib sees these macros.
cc=${CC-cc}
[[ -n $cc ]] || { echo "CC is empty; a compiler is required (see step 5)" >&2; exit 1; }

# Mirrors the `metis` target in cmake/EplinalgMumps.cmake. libmetis/ and
# GKlib/src/ both contain graph.c, so they must land in separate object
# directories or the census silently loses whichever gets clobbered.
flags=(-c -std=c99 -fPIC -w -fno-strict-aliasing
       -DIDXTYPEWIDTH=32 -DREALTYPEWIDTH=32 -DNDEBUG -DNDEBUG2
       -D_FILE_OFFSET_BITS=64 -DLINUX -DHAVE_EXECINFO_H -DHAVE_GETLINE
       -I"$dest/include" -I"$dest/GKlib/include" -I"$dest/libmetis")

# Compile the whole tree into $1, write every defined global to $2.
# Called plainly, never inside $(...), so a failed compile trips set -e
# instead of being swallowed by the command substitution's exit status.
census() {
    local obj=$1 out=$2 f
    rm -rf "$obj"; mkdir -p "$obj/gk" "$obj/lm"
    for f in "$dest"/GKlib/src/*.c; do "$cc" "${flags[@]}" "$f" -o "$obj/gk/$(basename "${f%.c}").o"; done
    for f in "$dest"/libmetis/*.c;   do "$cc" "${flags[@]}" "$f" -o "$obj/lm/$(basename "${f%.c}").o"; done
    nm --defined-only "$obj"/gk/*.o "$obj"/lm/*.o \
        | awk 'NF == 3 && $2 ~ /^[TDBRGWV]$/ { print $3 }' | sort -u > "$out"
}

say "enumerating GKlib exports with $cc"
census "$work/obj1" "$work/syms1.txt"
mapfile -t gk_exports < <(grep -E '^gk_' "$work/syms1.txt" || true)
[[ ${#gk_exports[@]} -gt 0 ]] || { echo "no gk_ exports found — GKlib layout changed?" >&2; exit 1; }
say "${#gk_exports[@]} gk_ exports + ${#GK_STRAYS[@]} unprefixed"

{
    printf '\n/* Private namespacing, matching libmetis/rename.h: every GKlib\n'
    printf '   global gets MUMPS_ prepended, so this copy cannot clash with\n'
    printf '   another GKlib in the link — gk_ is GKlib'"'"'s prefix, not ours.\n'
    printf '   Generated — see scripts/vendor_metis.sh. */\n'
    printf '%s\n' "${gk_exports[@]}" "${GK_STRAYS[@]}" | sort -u \
        | while read -r n; do printf '#define %-40s MUMPS_%s\n' "$n" "$n"; done
    printf '\n'
} > "$work/gk_rename.h"
awk -v block="$(cat "$work/gk_rename.h")" '
    /^#include <gk_types.h>/ && !done { print block; done = 1 }
    { print }
' "$dest/GKlib/include/GKlib.h" > "$work/GKlib.h.new"
grep -q '^#define gk_malloc  *MUMPS_gk_malloc$' "$work/GKlib.h.new" \
    || { echo "GKlib rename block not injected" >&2; exit 1; }
mv "$work/GKlib.h.new" "$dest/GKlib/include/GKlib.h"

# ── 6. verify ───────────────────────────────────────────────────────
# Every defined global must now carry one of the four private prefixes.
# A brand-new unprefixed name aborts and names itself.
say "verifying"
census "$work/obj2" "$work/syms2.txt"
bare=$(grep -vE '^(METIS_MUMPS_|metis_mumps_|libmetis_MUMPS_|MUMPS_)' "$work/syms2.txt" || true)
if [[ -n $bare ]]; then
    echo "ERROR: unprefixed globals survive the rename:" >&2
    printf '  %s\n' $bare >&2
    echo "Add each to LIBMETIS_STRAYS or GK_STRAYS in $0 and re-run." >&2
    exit 1
fi
say "$(wc -l < "$work/syms2.txt") globals, all private; 0 unprefixed"
say "done: $dest"
