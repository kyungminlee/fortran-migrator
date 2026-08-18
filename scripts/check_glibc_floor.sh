#!/bin/bash
#
# check_glibc_floor.sh — fail if the built archives cannot link against
# the glibc this script is running on.
#
# Usage:
#   scripts/check_glibc_floor.sh <dir>
#
#   <dir>  directory searched recursively for *.a
#
# Why this exists
# ---------------
# Everything this project ships is a *static* archive, so no glibc
# version is baked in at build time: a .a records only unversioned
# undefined symbol names, and the binding to a particular glibc version
# happens when the consumer links. That makes the release nearly
# host-agnostic — with one hole. A symbol the build host's headers chose
# but the consumer's glibc never exported at all cannot be bound at any
# version, and the consumer's link simply fails.
#
# Two such symbols shipped in v0.22.0, both from headers rather than
# source: __isoc23_sscanf (glibc >= 2.38, pulled in by _GNU_SOURCE in
# GKlib's io.c) and stat64 (glibc >= 2.33, pulled in by
# _FILE_OFFSET_BITS=64). Both are fixed; this script is what keeps the
# next one from shipping.
#
# So run it *on the oldest glibc you intend to support* — in CI that is a
# manylinux_2_28 container — and it answers exactly one question: is
# every symbol these archives reference either defined by the archives
# themselves, exported by this glibc, or supplied by a runtime the
# consumer brings anyway (libgfortran, libstdc++, MPI, OpenMP)?
#
set -euo pipefail
export LC_ALL=C

dir=${1:?usage: check_glibc_floor.sh <dir>}
[[ -d $dir ]] || { echo "not a directory: $dir" >&2; exit 2; }

work=$(mktemp -d); trap 'rm -rf "$work"' EXIT

mapfile -t archives < <(find "$dir" -name '*.a' | sort)
(( ${#archives[@]} > 0 )) || { echo "no .a files under $dir" >&2; exit 2; }

echo "glibc floor check"
echo "  glibc:    $(getconf GNU_LIBC_VERSION 2>/dev/null || ldd --version | head -1)"
echo "  archives: ${#archives[@]}"

# ── what the archives themselves provide / want ─────────────────────
nm --defined-only "${archives[@]}" 2>/dev/null \
    | awk 'NF == 3 && $2 ~ /^[A-Za-z]$/ { print $3 }' | sort -u > "$work/defined"
nm -u "${archives[@]}" 2>/dev/null \
    | awk 'NF == 2 && $1 == "U" { print $2 }' | sort -u > "$work/undefined"

# ── what the platform provides ──────────────────────────────────────
# glibc proper, plus the runtimes a consumer of this project links
# anyway. Absent ones are skipped: the allowlist below is the backstop,
# so a container without libgfortran still gives a usable answer.
: > "$work/provided"
# ldconfig's output is captured once rather than piped per soname: an awk
# that exits early would SIGPIPE the producer, and set -o pipefail turns
# that into a spurious failure.
ldconfig -p 2>/dev/null > "$work/ldcache" || : > "$work/ldcache"
for soname in libc.so.6 libm.so.6 libpthread.so.0 libdl.so.2 librt.so.1 \
              libgcc_s.so.1 libstdc++.so.6 libgfortran.so.5 libquadmath.so.0; do
    path=$(awk -v s="$soname" '$1 == s && /x86-64/ { print $NF; exit }' "$work/ldcache")
    [[ -n ${path:-} && -e $path ]] || continue
    echo "  runtime:  $path"
    nm -D --defined-only "$path" 2>/dev/null \
        | awk 'NF == 3 { sub(/@.*/, "", $3); print $3 }' >> "$work/provided"
done
# -lc is a linker script, not just libc.so.6: it pulls libc_nonshared.a
# alongside it. That archive is where glibc keeps the symbols it never
# exported dynamically, and its contents shrink over time — stat, fstat,
# lstat and mknod live there through 2.32 and only become real exports
# in 2.33. Reading libc.so.6 alone therefore reports a plain stat() as
# missing on exactly the old glibc we are trying to support, which is
# backwards. The static halves are counted too.
for stub in libc_nonshared.a libpthread_nonshared.a; do
    path=$( (gcc -print-file-name="$stub") 2>/dev/null || true)
    [[ -n ${path:-} && $path != "$stub" && -e $path ]] || continue
    echo "  runtime:  $path"
    nm --defined-only "$path" 2>/dev/null \
        | awk 'NF == 3 && $2 ~ /^[A-Za-z]$/ { print $3 }' >> "$work/provided"
done

sort -u -o "$work/provided" "$work/provided"

# ── symbols that are somebody else's problem ────────────────────────
# Not glibc's, so out of scope for a *glibc* floor: the consumer supplies
# these by linking a Fortran/C++/MPI/OpenMP runtime of their own choosing.
# blas_*_x_ are XBLAS's extended-precision entry points, and the last two
# come from the dynamic loader.
allow='^_gfortran|^_Z|^__gxx_|^__cxa_|^(MPI|PMPI|I_MPI)_|^(ompi|mpi|mpich)_|^(GOMP_|omp_|__kmpc)|^blas_[a-z0-9]+_x_$|^_GLOBAL_OFFSET_TABLE_$|^__tls_get_addr$'

comm -23 "$work/undefined" "$work/defined" \
    | comm -23 - "$work/provided" \
    | grep -Ev "$allow" > "$work/missing" || true

echo "  symbols:  $(wc -l < "$work/undefined") undefined, $(wc -l < "$work/defined") defined here"

if [[ -s $work/missing ]]; then
    echo
    echo "ERROR: $(wc -l < "$work/missing") symbol(s) this glibc does not export:" >&2
    while read -r sym; do
        printf '  %-28s referenced by %s\n' "$sym" \
            "$(nm -A -u "${archives[@]}" 2>/dev/null | awk -v s="$sym" '$NF == s { sub(/.*\//, "", $1); sub(/:$/, "", $1); print $1 }' | sort -u | paste -sd, - | cut -c1-90)" >&2
    done < "$work/missing"
    echo >&2
    echo "Each is either a header-driven symbol that needs pinning at the" >&2
    echo "build (see scripts/vendor_metis.sh and cmake/EplinalgMumps.cmake)" >&2
    echo "or a genuinely new dependency — if the latter, raise the floor" >&2
    echo "deliberately rather than by accident." >&2
    exit 1
fi

echo "  result:   OK — every reference resolves against this glibc"
