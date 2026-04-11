#!/usr/bin/env bash
# Compile sweep over a migrated multifloats BLACS source tree.
#
# Usage:  scripts/compile_blacs.sh [output-dir]
#
# Pipeline:
#   1. Build libmfc (the plain-C companion to multifloats: struct types
#      and runtime MPI handle registration).
#   2. Run the migrator over recipes/blacs.yaml with --target multifloats.
#   3. Compile every .c file in the output dir with mpicc.
#   4. Report ok/fail counts and write failing filenames to
#      /tmp/blacs_failures.txt.
#
# The script does NOT link a libddblacs.a archive -- that is the next
# step once compilation is clean. The intent here is to drive the
# c_migrator end-to-end and verify the migration produces buildable
# C source.

set -uo pipefail

REPO="$(cd "$(dirname "$0")"/.. && pwd)"
OUT_DIR="${1:-/tmp/mf_blacs}"
WORKERS="${WORKERS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

MFC_INC="$REPO/external/multifloats/include"
MFC_SRC="$REPO/external/multifloats/src"
MFC_OBJ_DIR="/tmp/mfc_build"
mkdir -p "$MFC_OBJ_DIR"

MPICC="${MPICC:-mpicc}"

# Drive the MPI wrapper through gcc/gfortran by default. The companion
# Fortran libraries (libblas, liblapack, libmultifloats) are built with
# gfortran via compile_lapack.sh, and using gcc-backed mpicc here keeps
# C ABI / name-mangling consistent with that toolchain. Override by
# setting OMPI_CC / OMPI_FC in the environment before invoking.
export OMPI_CC="${OMPI_CC:-gcc-15}"
export OMPI_FC="${OMPI_FC:-gfortran-15}"

echo "[1/4] Building libmfc..."
$MPICC -O2 -Wall -fPIC -std=c99 \
    -I "$MFC_INC" \
    -c "$MFC_SRC/multifloats_mpi.c" \
    -o "$MFC_OBJ_DIR/multifloats_mpi.o" || {
    echo "ERROR: failed to build libmfc"
    exit 1
}

echo "[2/4] Running migrator (recipes/blacs.yaml --target multifloats)..."
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"
# Use src-multifloats/pyengine as the canonical pyengine for multifloats
# work; the root-level pyengine/ is the older KIND-only version.
(cd "$REPO/src-multifloats" && uv run --project .. python -m pyengine migrate \
    "$REPO/recipes/blacs.yaml" "$OUT_DIR" \
    --target multifloats) || {
    echo "ERROR: migration failed"
    exit 1
}

echo "[3/4] Compiling migrated BLACS sources with $MPICC..."
mkdir -p "$OUT_DIR/objs"
: > /tmp/blacs_failures.txt
: > /tmp/blacs_compile_results.txt

# BLACS compile flags (mirrors stock ScaLAPACK BLACS Bmake.MPI):
#   -DCSAMEF77=1   integer Fortran types match C int
#   -DInt=int      Int typedef
#   -DUseMpi2=1    enable MPI-2 features
#   -DAdd_         Fortran name mangling (lowercase + trailing underscore)
BLACS_DEFS="-DCSAMEF77=1 -DInt=int -DUseMpi2=1 -DAdd_"

compile_one() {
    f="$1"
    out_dir="$2"
    obj="$out_dir/objs/$(basename "$f").o"
    err=$($MPICC -O2 -Wall -std=c99 \
        -I "$MFC_INC" \
        -I "$out_dir" \
        $BLACS_DEFS \
        -c "$f" -o "$obj" 2>&1)
    rc=$?
    if [ $rc -eq 0 ]; then
        echo "OK $(basename "$f")"
    else
        echo "FAIL $(basename "$f")"
        echo "=== $(basename "$f") ===" >> /tmp/blacs_failures.txt
        echo "$err" >> /tmp/blacs_failures.txt
    fi
}
export -f compile_one
export MPICC MFC_INC BLACS_DEFS

find "$OUT_DIR" -maxdepth 1 -name '*.c' -print0 \
    | xargs -0 -n1 -P "$WORKERS" -I {} \
        bash -c 'compile_one "$@"' _ {} "$OUT_DIR" \
    > /tmp/blacs_compile_results.txt

ok=$(grep -c '^OK ' /tmp/blacs_compile_results.txt || echo 0)
fail=$(grep -c '^FAIL ' /tmp/blacs_compile_results.txt || echo 0)

echo "[4/4] Done. $ok ok / $fail failed"
echo "    full log:    /tmp/blacs_compile_results.txt"
echo "    failures:    /tmp/blacs_failures.txt"
