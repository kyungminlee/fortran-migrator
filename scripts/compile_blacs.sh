#!/usr/bin/env bash
# Compile sweep over a migrated multifloats BLACS source tree.
#
# Usage:  scripts/compile_blacs.sh [output-dir]
#
# Pipeline:
#   1. Compile multifloats_mpi.cpp (MPI handle registration).
#   2. Run the migrator over recipes/blacs.yaml with --target multifloats.
#   3. Compile every .c file as C++ with mpicxx + the bridge header.
#   4. Report ok/fail counts.

set -uo pipefail

REPO="$(cd "$(dirname "$0")"/.. && pwd)"
OUT_DIR="${1:-/tmp/mf_blacs}"
WORKERS="${WORKERS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

MF_BRIDGE_INC="$REPO/external/multifloats/include"
MF_HH_DIR="$REPO/external"
MF_MPI_SRC="$REPO/external/multifloats/src/multifloats_mpi.cpp"

MPICXX="${MPICXX:-mpicxx}"

export OMPI_CC="${OMPI_CC:-gcc-15}"
export OMPI_CXX="${OMPI_CXX:-g++-15}"
export OMPI_FC="${OMPI_FC:-gfortran-15}"

echo "[1/4] Compiling multifloats_mpi.cpp..."
$MPICXX -O2 -std=c++17 \
    -I "$MF_BRIDGE_INC" -I "$MF_HH_DIR" \
    -c "$MF_MPI_SRC" \
    -o "/tmp/multifloats_mpi.o" || {
    echo "ERROR: failed to compile multifloats_mpi.cpp"
    exit 1
}

echo "[2/4] Running migrator (recipes/blacs.yaml --target multifloats)..."
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"
(cd "$REPO/src" && uv run --project .. python -m migrator migrate \
    "$REPO/recipes/blacs.yaml" "$OUT_DIR" \
    --target multifloats) || {
    echo "ERROR: migration failed"
    exit 1
}

# Copy the bridge header into the output directory so overrides can find it
cp "$MF_BRIDGE_INC/multifloats_bridge.h" "$OUT_DIR/"

echo "[3/4] Compiling migrated BLACS sources as C++ with $MPICXX..."
mkdir -p "$OUT_DIR/objs"

BLACS_DEFS="-DCSAMEF77=1 -DInt=int -DUseMpi2=1 -DAdd_"

compile_one() {
    f="$1"
    out_dir="$2"
    obj="$out_dir/objs/$(basename "$f").o"
    err=$($MPICXX -O2 -std=c++17 \
        -fpermissive -Wno-write-strings \
        -I "$MF_BRIDGE_INC" \
        -I "$MF_HH_DIR" \
        -I "$out_dir" \
        -include multifloats_bridge.h \
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
export MPICXX MF_BRIDGE_INC MF_HH_DIR BLACS_DEFS

: > /tmp/blacs_failures.txt

find "$OUT_DIR" -maxdepth 1 -name '*.c' -print0 \
    | xargs -0 -n1 -P "$WORKERS" -I {} \
        bash -c 'compile_one "$@"' _ {} "$OUT_DIR" \
    > /tmp/blacs_compile_results.txt

ok=$(grep -c '^OK ' /tmp/blacs_compile_results.txt || echo 0)
fail=$(grep -c '^FAIL ' /tmp/blacs_compile_results.txt || echo 0)

echo "[4/4] Done. $ok ok / $fail failed"
echo "    full log:    /tmp/blacs_compile_results.txt"
echo "    failures:    /tmp/blacs_failures.txt"
