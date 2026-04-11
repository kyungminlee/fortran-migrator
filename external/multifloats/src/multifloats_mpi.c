/* multifloats_mpi.c -- runtime registration of multifloats MPI handles.
 *
 * Defines the global MPI_Datatype and MPI_Op handles declared extern in
 * multifloats_c.h. multifloats_mpi_init() performs one-time registration
 * using standard MPI primitives; it is called from the migrated BLACS
 * grid-init entry points so user code does not need to call it.
 */
#include "multifloats_c.h"

MPI_Datatype MPI_FLOAT64X2    = MPI_DATATYPE_NULL;
MPI_Datatype MPI_COMPLEX128X2 = MPI_DATATYPE_NULL;

MPI_Op MPI_DD_SUM = MPI_OP_NULL;
MPI_Op MPI_ZZ_SUM = MPI_OP_NULL;
MPI_Op MPI_DD_AMX = MPI_OP_NULL;
MPI_Op MPI_DD_AMN = MPI_OP_NULL;
MPI_Op MPI_ZZ_AMX = MPI_OP_NULL;
MPI_Op MPI_ZZ_AMN = MPI_OP_NULL;

/* --------------------------------------------------------------------- */
/* User-op callbacks                                                      */
/* --------------------------------------------------------------------- */

static void mf_dd_sum_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    (void)dt;
    float64x2_t *a = (float64x2_t *)in;
    float64x2_t *b = (float64x2_t *)inout;
    for (int i = 0; i < *len; ++i) {
        b[i] = mf_add(a[i], b[i]);
    }
}

static void mf_zz_sum_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    (void)dt;
    complex128x2_t *a = (complex128x2_t *)in;
    complex128x2_t *b = (complex128x2_t *)inout;
    for (int i = 0; i < *len; ++i) {
        b[i] = mf_cadd(a[i], b[i]);
    }
}

static void mf_dd_amx_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    (void)dt;
    float64x2_t *a = (float64x2_t *)in;
    float64x2_t *b = (float64x2_t *)inout;
    for (int i = 0; i < *len; ++i) {
        if (mf_cmp(mf_abs(a[i]), mf_abs(b[i])) > 0) {
            b[i] = a[i];
        }
    }
}

static void mf_dd_amn_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    (void)dt;
    float64x2_t *a = (float64x2_t *)in;
    float64x2_t *b = (float64x2_t *)inout;
    for (int i = 0; i < *len; ++i) {
        if (mf_cmp(mf_abs(a[i]), mf_abs(b[i])) < 0) {
            b[i] = a[i];
        }
    }
}

static void mf_zz_amx_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    (void)dt;
    complex128x2_t *a = (complex128x2_t *)in;
    complex128x2_t *b = (complex128x2_t *)inout;
    for (int i = 0; i < *len; ++i) {
        if (mf_cmp(mf_cabs1(a[i]), mf_cabs1(b[i])) > 0) {
            b[i] = a[i];
        }
    }
}

static void mf_zz_amn_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    (void)dt;
    complex128x2_t *a = (complex128x2_t *)in;
    complex128x2_t *b = (complex128x2_t *)inout;
    for (int i = 0; i < *len; ++i) {
        if (mf_cmp(mf_cabs1(a[i]), mf_cabs1(b[i])) < 0) {
            b[i] = a[i];
        }
    }
}

/* --------------------------------------------------------------------- */
/* One-time registration                                                  */
/* --------------------------------------------------------------------- */

void multifloats_mpi_init(void) {
    /* BLACS is not normally called from multiple threads concurrently and
     * the first BLACS call happens after MPI_Init has returned, so a plain
     * static flag is sufficient. Idempotence: the function is safe to call
     * any number of times and all subsequent calls after the first are
     * single-branch no-ops. */
    static int initialized = 0;
    if (initialized) {
        return;
    }

    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_FLOAT64X2);
    MPI_Type_commit(&MPI_FLOAT64X2);

    MPI_Type_contiguous(4, MPI_DOUBLE, &MPI_COMPLEX128X2);
    MPI_Type_commit(&MPI_COMPLEX128X2);

    /* commute=1: all six ops are commutative, matching stock MPI_SUM /
     * MPI_MAX / MPI_MIN semantics. Associativity is only approximate --
     * same as stock double -- but BLACS already warns about this. */
    MPI_Op_create(mf_dd_sum_fn, 1, &MPI_DD_SUM);
    MPI_Op_create(mf_zz_sum_fn, 1, &MPI_ZZ_SUM);
    MPI_Op_create(mf_dd_amx_fn, 1, &MPI_DD_AMX);
    MPI_Op_create(mf_dd_amn_fn, 1, &MPI_DD_AMN);
    MPI_Op_create(mf_zz_amx_fn, 1, &MPI_ZZ_AMX);
    MPI_Op_create(mf_zz_amn_fn, 1, &MPI_ZZ_AMN);

    initialized = 1;
}
