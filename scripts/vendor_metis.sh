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
#   CC       compiler used for the verification build (default: cc).
#            Set CC=  (empty) to skip verification entirely.
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
#   3. Internals upstream forgot to list in rename.h (LIBMETIS_STRAYS)
#      and GKlib exports that escaped GKlib's own gk_ convention
#      (GK_STRAYS). These are appended as extra #defines rather than
#      edited in place, matching how upstream handles the rest.
#
#   4. metis.h's IDXTYPEWIDTH / REALTYPEWIDTH self-defaults, which 5.2.1
#      comments out, are restored to 32 — the width this project builds
#      and the 5.1.0 behavior. Without them every translation unit that
#      includes the *installed* metis.h would have to supply -D itself.
#
# Everything else is upstream verbatim; `diff -r` against a pristine
# checkout with the rename reversed comes back empty.
#
# Verification (unless CC is empty): the result is compiled with this
# project's flags and every defined global in the resulting objects is
# checked to carry one of the four private prefixes. A brand-new
# unprefixed name aborts the run and names itself — add it to
# LIBMETIS_STRAYS or GK_STRAYS and re-run.
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

# ── 3. restore the width self-defaults ──────────────────────────────
sed -i -e 's|^//#define IDXTYPEWIDTH 32$|#define IDXTYPEWIDTH 32|' \
       -e 's|^//#define REALTYPEWIDTH 32$|#define REALTYPEWIDTH 32|' \
       "$dest/include/metis.h"
grep -q '^#define IDXTYPEWIDTH 32$'  "$dest/include/metis.h" || { echo "IDXTYPEWIDTH patch failed" >&2; exit 1; }
grep -q '^#define REALTYPEWIDTH 32$' "$dest/include/metis.h" || { echo "REALTYPEWIDTH patch failed" >&2; exit 1; }

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

# GKlib's own exports carry gk_; these do not. GKlib.h is the one header
# every GKlib and libmetis translation unit sees, and the block goes
# immediately before <gk_proto.h> so it covers the declarations too.
{
    printf '\n/* Private namespacing, matching libmetis/rename.h: these GKlib\n'
    printf '   exports lack GKlib'"'"'s own gk_ prefix and would otherwise be\n'
    printf '   bare globals in the archive. See scripts/vendor_metis.sh. */\n'
    for n in "${GK_STRAYS[@]}"; do
        printf '#define %-32s gk_MUMPS_%s\n' "$n" "$n"
    done
    printf '\n'
} > "$work/gk_strays.h"
awk -v strays="$(cat "$work/gk_strays.h")" '
    /^#include <gk_proto.h>/ && !done { print strays; done = 1 }
    { print }
' "$dest/GKlib/include/GKlib.h" > "$work/GKlib.h.new"
mv "$work/GKlib.h.new" "$dest/GKlib/include/GKlib.h"

# gkregex.c defines its no-op stub *above* its own #include "GKlib.h",
# so the macro above cannot reach it; rename the definition in place.
sed -i 's/^void gkfooo()/void gk_MUMPS_gkfooo()/' "$dest/GKlib/src/gkregex.c"
grep -q '^void gk_MUMPS_gkfooo()' "$dest/GKlib/src/gkregex.c" || { echo "gkfooo patch failed" >&2; exit 1; }

# ── 5. verify ───────────────────────────────────────────────────────
cc=${CC-cc}
if [[ -z $cc ]]; then
    say "CC empty — skipping verification build"
    say "done: $dest"
    exit 0
fi

say "verifying with $cc"
# Mirrors the `metis` target in cmake/EplinalgMumps.cmake. libmetis/ and
# GKlib/src/ both contain graph.c, so they must land in separate object
# directories or the census silently loses whatever gets clobbered.
flags=(-c -std=c99 -fPIC -w -fno-strict-aliasing
       -DIDXTYPEWIDTH=32 -DREALTYPEWIDTH=32 -DNDEBUG -DNDEBUG2
       -D_FILE_OFFSET_BITS=64 -DLINUX -DHAVE_EXECINFO_H -DHAVE_GETLINE
       -I"$dest/include" -I"$dest/GKlib/include" -I"$dest/libmetis")
obj=$work/obj
rm -rf "$obj"; mkdir -p "$obj/gk" "$obj/lm"
for f in "$dest"/GKlib/src/*.c; do "$cc" "${flags[@]}" "$f" -o "$obj/gk/$(basename "${f%.c}").o"; done
for f in "$dest"/libmetis/*.c;   do "$cc" "${flags[@]}" "$f" -o "$obj/lm/$(basename "${f%.c}").o"; done

bare=$(nm --defined-only "$obj"/gk/*.o "$obj"/lm/*.o \
       | awk 'NF == 3 && $2 ~ /^[TDBRGWV]$/ { print $3 }' \
       | grep -vE '^(METIS_MUMPS_|metis_mumps_|libmetis_MUMPS_|gk_)' | sort -u || true)
if [[ -n $bare ]]; then
    echo "ERROR: unprefixed globals survive the rename:" >&2
    printf '  %s\n' $bare >&2
    echo "Add each to LIBMETIS_STRAYS or GK_STRAYS in $0 and re-run." >&2
    exit 1
fi
say "$(ls "$obj"/gk/*.o "$obj"/lm/*.o | wc -l) objects, 0 unprefixed globals"
say "done: $dest"
