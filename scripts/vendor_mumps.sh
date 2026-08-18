#!/bin/bash
#
# vendor_mumps.sh — regenerate extern/MUMPS_5.9.1 from upstream.
#
# The vendored MUMPS tree is upstream's tarball verbatim except for one
# mechanical edit: the eight METIS entry points MUMPS calls are rewritten
# to the private names the vendored METIS actually exports. METIS is
# namespaced (see scripts/vendor_metis.sh) so the ordering archive cannot
# clash with a system METIS the final application might also link; MUMPS
# is the only consumer, so it has to be taught the private names.
#
# The tree is checked into git already renamed. This script is how it is
# produced, so the edit stays mechanical, auditable, and repeatable when
# upstream moves. It is a maintainer tool — the build never runs it.
#
# Usage:
#   scripts/vendor_mumps.sh [<dest>]
#
#   <dest>  output directory (default: extern/MUMPS_5.9.1). Removed and
#           recreated. Nothing else in the repo is touched — re-pointing
#           codegen/recipes/mumps.yaml, codegen/migrator/staging.py and
#           the docs at a new version is a manual follow-up.
#
# Environment:
#   WORKDIR  scratch directory for the download and the verification
#            copy (default: a mktemp dir, removed on exit). Point it at
#            a persistent path to avoid re-downloading while iterating.
#   TARBALL  pre-downloaded MUMPS_<version>.tar.gz to use instead of
#            fetching. Still checksummed.
#
# The rename, in full:
#
#   METIS_NodeND            -> METIS_MUMPS_NodeND
#   METIS_NODEND            -> METIS_MUMPS_NODEND
#   METIS_NODEWND           -> METIS_MUMPS_NODEWND
#   METIS_PartGraphKway     -> METIS_MUMPS_PartGraphKway
#   METIS_SetDefaultOptions -> METIS_MUMPS_SetDefaultOptions
#   METIS_SETDEFAULTOPTIONS -> METIS_MUMPS_SETDEFAULTOPTIONS
#   metis_nodend            -> metis_mumps_nodend
#   metis_setdefaultoptions -> metis_mumps_setdefaultoptions
#
# Only those eight, and only whole-token. MUMPS's own identifiers share
# the METIS_ prefix without being METIS symbols — METIS_OPTIONS,
# METIS_IDX_SIZE, METIS_NOPTIONS, METIS_LIB, metis_lib_name, and the
# mixed-width wrappers METIS_NODEND_MIXEDto32 / METIS_NODEWND_MIXEDto64 —
# and a blanket prefix rewrite would break every one of them. The four
# stems below are matched case-insensitively and each spelling found is
# rewritten in its own casing, so a Fortran wrapper that starts emitting
# a new casing is picked up automatically.
#
# METIS_NODEND / METIS_NODEWND are the METIS 4 entry points; MUMPS still
# declares them behind its metis4 guards. They are renamed for symmetry
# even though the vendored METIS 5 does not define them.
#
# Verification: the rename is undone on a scratch copy with the two
# blanket rules METIS_MUMPS_ -> METIS_ and metis_mumps_ -> metis_, and
# the result must be byte-identical to the extracted tarball. That is a
# proof of the complete delta — any accidental edit anywhere in the tree
# shows up as a diff. The occurrence count is checked too, so a rule that
# silently matches nothing cannot pass.
#
set -euo pipefail

MUMPS_VERSION=5.9.1
MUMPS_URL=https://mumps-solver.org/MUMPS_${MUMPS_VERSION}.tar.gz
MUMPS_SHA256=659c9b57646b5a003ac618baa1faf9dd2044e46c732b3daaccbc7158003e1b46

# METIS entry points MUMPS calls. Matched case-insensitively; every
# casing present in the tarball is renamed in place.
METIS_STEMS=(NodeND NodeWND PartGraphKway SetDefaultOptions)

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
dest=${1:-$root/extern/MUMPS_${MUMPS_VERSION}}
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
tar=${TARBALL:-$work/MUMPS_${MUMPS_VERSION}.tar.gz}
if [[ ! -f $tar ]]; then
    say "downloading $MUMPS_URL"
    curl -fsSL "$MUMPS_URL" -o "$tar"
fi
have=$(sha256sum "$tar" | cut -d' ' -f1)
[[ $have == "$MUMPS_SHA256" ]] || {
    echo "$tar sha256 $have, expected $MUMPS_SHA256" >&2; exit 1; }
say "sha256 ok"

pristine=$work/pristine
rm -rf "$pristine"; mkdir -p "$pristine"
tar -xzf "$tar" -C "$pristine"
[[ -d $pristine/MUMPS_${MUMPS_VERSION} ]] || {
    echo "tarball does not contain MUMPS_${MUMPS_VERSION}/" >&2; exit 1; }
pristine=$pristine/MUMPS_${MUMPS_VERSION}

# The prefix must not already be upstream's; if it ever is, the reverse
# rules below stop being a faithful inverse and the proof is worthless.
if grep -rqI 'METIS_MUMPS_\|metis_mumps_' "$pristine"; then
    echo "upstream already uses the METIS_MUMPS_ prefix — revisit this script" >&2
    exit 1
fi

# ── 2. copy the tree verbatim ───────────────────────────────────────
# Unpruned on purpose: the migrator reads src/ and include/, the build
# takes libseq/ and PORD/, and the rest (doc, examples, Make.inc) is
# provenance that costs little and answers questions later.
say "assembling $dest"
rm -rf "$dest"
cp -a "$pristine" "$dest"

# ── 3. rename ───────────────────────────────────────────────────────
sedprog=$work/rename.sed
: > "$sedprog"
for stem in "${METIS_STEMS[@]}"; do
    # Every distinct casing of METIS_<stem> that actually occurs.
    ci=$(printf '%s' "$stem" | sed 's/./[\U&\L&]/g')
    mapfile -t spellings < <(
        grep -rhoIE "\\b[Mm][Ee][Tt][Ii][Ss]_$ci\\b" "$dest" | sort -u)
    [[ ${#spellings[@]} -gt 0 ]] || {
        echo "no occurrence of METIS_$stem — stem list is stale" >&2; exit 1; }
    for name in "${spellings[@]}"; do
        case $name in
            METIS_*) new=METIS_MUMPS_${name#METIS_} ;;
            metis_*) new=metis_mumps_${name#metis_} ;;
            *) echo "unhandled casing '$name' — extend the case above" >&2; exit 1 ;;
        esac
        printf 's/\\b%s\\b/%s/g\n' "$name" "$new" >> "$sedprog"
        say "$name -> $new"
    done
done

# Only the files that actually match, and never a binary: the tree ships
# a PDF user guide and a MATLAB .mat, and a blanket sed over those is a
# corruption waiting to happen that the round-trip below cannot see.
pattern=$(sed -n 's|^s/\\b\([^\\]*\)\\b.*|\1|p' "$sedprog" | paste -sd'|')
mapfile -t targets < <(grep -rlIE "\\b($pattern)\\b" "$dest")
say "${#targets[@]} files touched"
before=$(grep -rhoIE "\\b($pattern)\\b" "$dest" | wc -l)
sed -i -f "$sedprog" "${targets[@]}"

# ── 4. verify ───────────────────────────────────────────────────────
left=$(grep -rhoIE "\\b($pattern)\\b" "$dest" | sort | uniq -c || true)
[[ -z $left ]] || { echo "ERROR: unrenamed occurrences survive:" >&2
                    printf '%s\n' "$left" >&2; exit 1; }
after=$(grep -rhoI 'METIS_MUMPS_\|metis_mumps_' "$dest" | wc -l)
[[ $before -eq $after ]] || {
    echo "ERROR: renamed $after occurrences, found $before before the edit" >&2
    exit 1; }

rev=$work/reverse
rm -rf "$rev"; cp -a "$dest" "$rev"
grep -rlIE 'METIS_MUMPS_|metis_mumps_' "$rev" \
    | xargs sed -i -e 's/METIS_MUMPS_/METIS_/g' -e 's/metis_mumps_/metis_/g'
if ! diff -rq "$pristine" "$rev" > "$work/diff.txt" 2>&1; then
    echo "ERROR: the tree is not upstream-plus-the-rename:" >&2
    cat "$work/diff.txt" >&2
    exit 1
fi
say "$before occurrences renamed; reverse round-trip is byte-identical to upstream"
say "done: $dest"
