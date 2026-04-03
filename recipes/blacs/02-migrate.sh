#!/usr/bin/env bash
# 02-migrate.sh — Clone BLACS type-specific C files for extended precision.
#
# Creates new q* (extended real) and x* (extended complex) variants
# by cloning d* and z* files with mechanical text substitution.
# Also copies infrastructure files and generates a modified Bdef.h.
#
# Usage: ./recipes/blacs/02-migrate.sh [--kind 10|16]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00-setup.sh" "$@"

echo "=== Migrating BLACS to KIND=$TARGET_KIND ==="

RP="$REAL_PREFIX"
CP="$COMPLEX_PREFIX"
RPU="$REAL_PREFIX_UPPER"
CPU="$COMPLEX_PREFIX_UPPER"

# Clean output
rm -rf "${OUTPUT_DIR:?}"/*

cd "$BLACS_SRC"

# --------------------------------------------------------------------------
# Step 1: Copy all original files
# --------------------------------------------------------------------------
echo "Copying original files..."
cp -a *.c *.h "$OUTPUT_DIR/" 2>/dev/null || true

# --------------------------------------------------------------------------
# Step 2: Clone d-variant → q/e-variant (real extended)
# --------------------------------------------------------------------------
echo "Cloning d-variants → ${RP}-variants..."

clone_d_to_ext() {
    local src="$1"
    local dst
    # Replace leading 'd' in filename, or 'BI_d' prefix
    dst=$(echo "$src" | sed "s/^d/${RP}/" | sed "s/^BI_d/BI_${RP}/")

    # Use -E for extended regex. macOS sed doesn't support \b, so we
    # match word boundaries with [^a-zA-Z_] context and run twice to
    # catch adjacent matches.
    sed -E \
        -e "s/(^|[^a-zA-Z_])double([^a-zA-Z_]|$)/\1QREAL\2/g" \
        -e "s/(^|[^a-zA-Z_])double([^a-zA-Z_]|$)/\1QREAL\2/g" \
        -e "s/MPI_DOUBLE/MPI_QREAL/g" \
        -e "s/Cd([a-z])/C${RP}\1/g" \
        -e "s/BI_d([a-z])/BI_${RP}\1/g" \
        -e "s/([^a-zA-Z_])d(ge[srvb][sdvr]2d)/\1${RP}\2/g" \
        -e "s/([^a-zA-Z_])d(tr[srvb][sdvr]2d)/\1${RP}\2/g" \
        -e "s/([^a-zA-Z_])d(g[sa][umx][mnx]2d)/\1${RP}\2/g" \
        -e "s/([^a-zA-Z_])d(gsum2d)/\1${RP}\2/g" \
        -e "s/^d(ge[srvb][sdvr]2d)/${RP}\1/g" \
        -e "s/^d(tr[srvb][sdvr]2d)/${RP}\1/g" \
        -e "s/^d(g[sa][umx][mnx]2d)/${RP}\1/g" \
        -e "s/^d(gsum2d)/${RP}\1/g" \
        -e "s/\"double precision\"/\"extended precision\"/g" \
        "$BLACS_SRC/$src" > "$OUTPUT_DIR/$dst"

    echo "  $src → $dst"
}

# User-facing d-routines
for f in d{ge,tr}{sd,rv,bs,br}2d_.c d{gsum,gamx,gamn}2d_.c; do
    [[ -f "$f" ]] && clone_d_to_ext "$f"
done

# BI_d vector ops
for f in BI_dvv{sum,amn,amn2,amx,amx2}.c; do
    [[ -f "$f" ]] && clone_d_to_ext "$f"
done

# BI_d MPI wrappers
for f in BI_dMPI_*.c; do
    [[ -f "$f" ]] && clone_d_to_ext "$f"
done

# BI_d matrix copy
for f in BI_dmvcopy.c BI_dvmcopy.c; do
    [[ -f "$f" ]] && clone_d_to_ext "$f"
done

# --------------------------------------------------------------------------
# Step 3: Clone z-variant → x/y-variant (complex extended)
# --------------------------------------------------------------------------
echo ""
echo "Cloning z-variants → ${CP}-variants..."

clone_z_to_ext() {
    local src="$1"
    local dst
    dst=$(echo "$src" | sed "s/^z/${CP}/" | sed "s/^BI_z/BI_${CP}/")

    sed -E \
        -e "s/(^|[^a-zA-Z_])DCOMPLEX([^a-zA-Z_]|$)/\1QCOMPLEX\2/g" \
        -e "s/(^|[^a-zA-Z_])DCOMPLEX([^a-zA-Z_]|$)/\1QCOMPLEX\2/g" \
        -e "s/(^|[^a-zA-Z_])SCOMPLEX([^a-zA-Z_]|$)/\1QCOMPLEX\2/g" \
        -e "s/(^|[^a-zA-Z_])SCOMPLEX([^a-zA-Z_]|$)/\1QCOMPLEX\2/g" \
        -e "s/(^|[^a-zA-Z_])double([^a-zA-Z_]|$)/\1QREAL\2/g" \
        -e "s/(^|[^a-zA-Z_])double([^a-zA-Z_]|$)/\1QREAL\2/g" \
        -e "s/MPI_DOUBLE_COMPLEX/MPI_QCOMPLEX/g" \
        -e "s/MPI_DOUBLE/MPI_QREAL/g" \
        -e "s/MPI_COMPLEX([^a-zA-Z_0-9])/MPI_QCOMPLEX\1/g" \
        -e "s/Cz([a-z])/C${CP}\1/g" \
        -e "s/BI_z([a-z])/BI_${CP}\1/g" \
        -e "s/BI_d([a-z])/BI_${RP}\1/g" \
        -e "s/([^a-zA-Z_])z(ge[srvb][sdvr]2d)/\1${CP}\2/g" \
        -e "s/([^a-zA-Z_])z(tr[srvb][sdvr]2d)/\1${CP}\2/g" \
        -e "s/([^a-zA-Z_])z(g[sa][umx][mnx]2d)/\1${CP}\2/g" \
        -e "s/([^a-zA-Z_])z(gsum2d)/\1${CP}\2/g" \
        -e "s/^z(ge[srvb][sdvr]2d)/${CP}\1/g" \
        -e "s/^z(tr[srvb][sdvr]2d)/${CP}\1/g" \
        -e "s/^z(g[sa][umx][mnx]2d)/${CP}\1/g" \
        -e "s/^z(gsum2d)/${CP}\1/g" \
        -e "s/\"double complex\"/\"extended complex\"/g" \
        "$BLACS_SRC/$src" > "$OUTPUT_DIR/$dst"

    echo "  $src → $dst"
}

# User-facing z-routines
for f in z{ge,tr}{sd,rv,bs,br}2d_.c z{gsum,gamx,gamn}2d_.c; do
    [[ -f "$f" ]] && clone_z_to_ext "$f"
done

# BI_z vector ops
for f in BI_zvv{sum,amn,amn2,amx,amx2}.c; do
    [[ -f "$f" ]] && clone_z_to_ext "$f"
done

# BI_z MPI wrappers
for f in BI_zMPI_*.c; do
    [[ -f "$f" ]] && clone_z_to_ext "$f"
done

# --------------------------------------------------------------------------
# Step 4: Generate extended Bdef.h
# --------------------------------------------------------------------------
echo ""
echo "Generating extended Bdef.h..."

# Add extended-precision type definitions before the #endif
cat >> "$OUTPUT_DIR/Bdef.h" << 'HEADER_EOF'

/*
 * ========================================================================
 *     EXTENDED PRECISION TYPE DEFINITIONS (auto-generated by migrator)
 * ========================================================================
 */
#ifdef HAVE_FLOAT128
typedef __float128 QREAL;
typedef struct {__float128 r, i;} QCOMPLEX;
#else
typedef long double QREAL;
typedef struct {long double r, i;} QCOMPLEX;
#endif

/* Extended-precision data type constants */
#define QUAD      8
#define COMPLEX32 9

/* Extended-precision copy macros */
void BI_qmvcopy(Int m, Int n, QREAL *A, Int lda, QREAL *buff);
void BI_qvmcopy(Int m, Int n, QREAL *A, Int lda, QREAL *buff);
#define BI_xmvcopy(m, n, A, lda, buff) \
        BI_qmvcopy(2*(m), (n), (QREAL *) (A), 2*(lda), (QREAL *) (buff))
#define BI_xvmcopy(m, n, A, lda, buff) \
        BI_qvmcopy(2*(m), (n), (QREAL *) (A), 2*(lda), (QREAL *) (buff))

/* Extended-precision MPI datatypes (must be initialized at runtime) */
extern MPI_Datatype MPI_QREAL;
extern MPI_Datatype MPI_QCOMPLEX;
void BI_InitQuadMpiTypes(void);

HEADER_EOF

echo "  Added extended-precision definitions to Bdef.h"

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
echo ""
orig_count=$(ls "$BLACS_SRC"/*.c 2>/dev/null | wc -l | tr -d ' ')
new_count=$(ls "$OUTPUT_DIR"/*.c 2>/dev/null | wc -l | tr -d ' ')
added=$((new_count - orig_count))

echo "=== Migration Summary ==="
echo "  Original files: $orig_count"
echo "  Output files:   $new_count (+$added new)"
echo "  Output dir:     $OUTPUT_DIR"
echo ""
echo "New files:"
diff <(ls "$BLACS_SRC"/*.c | xargs -I{} basename {} | sort) \
     <(ls "$OUTPUT_DIR"/*.c | xargs -I{} basename {} | sort) \
     | grep '^>' | sed 's/^> /  /' || true
