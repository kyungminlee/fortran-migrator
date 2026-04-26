/* multifloats_mpi.cpp -- runtime registration of multifloats MPI handles.
 *
 * Uses the C++ multifloats::float64x2 operators (+, abs via < 0 ? -x : x,
 * comparison) so no manual Knuth two-sum is needed.
 *
 * Compiled as C++; all exported symbols use extern "C" linkage so they
 * are callable from both C and Fortran (through the migrated BLACS
 * wrappers).
 */
#include "multifloats_bridge.h"

extern "C" {

MPI_Datatype MPI_FLOAT64X2    = MPI_DATATYPE_NULL;
MPI_Datatype MPI_COMPLEX128X2 = MPI_DATATYPE_NULL;

MPI_Op MPI_DD_SUM = MPI_OP_NULL;
MPI_Op MPI_ZZ_SUM = MPI_OP_NULL;
MPI_Op MPI_DD_AMX = MPI_OP_NULL;
MPI_Op MPI_DD_AMN = MPI_OP_NULL;
MPI_Op MPI_ZZ_AMX = MPI_OP_NULL;
MPI_Op MPI_ZZ_AMN = MPI_OP_NULL;

} /* extern "C" */

/* ---- User-op callbacks ------------------------------------------ */

static void dd_sum_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<float64x2 *>(in);
    auto *b = static_cast<float64x2 *>(inout);
    for (int i = 0; i < *len; ++i) b[i] = b[i] + a[i];
}

static void zz_sum_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<complex64x2 *>(in);
    auto *b = static_cast<complex64x2 *>(inout);
    for (int i = 0; i < *len; ++i) {
        b[i].re = b[i].re + a[i].re;
        b[i].im = b[i].im + a[i].im;
    }
}

static void dd_amx_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<float64x2 *>(in);
    auto *b = static_cast<float64x2 *>(inout);
    for (int i = 0; i < *len; ++i)
        if (mf_abs(a[i]) > mf_abs(b[i])) b[i] = a[i];
}

static void dd_amn_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<float64x2 *>(in);
    auto *b = static_cast<float64x2 *>(inout);
    for (int i = 0; i < *len; ++i)
        if (mf_abs(a[i]) < mf_abs(b[i])) b[i] = a[i];
}

static void zz_amx_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<complex64x2 *>(in);
    auto *b = static_cast<complex64x2 *>(inout);
    for (int i = 0; i < *len; ++i)
        if (mf_cabs1(a[i]) > mf_cabs1(b[i])) b[i] = a[i];
}

static void zz_amn_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<complex64x2 *>(in);
    auto *b = static_cast<complex64x2 *>(inout);
    for (int i = 0; i < *len; ++i)
        if (mf_cabs1(a[i]) < mf_cabs1(b[i])) b[i] = a[i];
}

/* ---- One-time registration -------------------------------------- */

extern "C" void multifloats_mpi_init(void) {
    static int initialized = 0;
    if (initialized) return;

    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_FLOAT64X2);
    MPI_Type_commit(&MPI_FLOAT64X2);
    MPI_Type_contiguous(4, MPI_DOUBLE, &MPI_COMPLEX128X2);
    MPI_Type_commit(&MPI_COMPLEX128X2);

    MPI_Op_create(dd_sum_fn, 1, &MPI_DD_SUM);
    MPI_Op_create(zz_sum_fn, 1, &MPI_ZZ_SUM);
    MPI_Op_create(dd_amx_fn, 1, &MPI_DD_AMX);
    MPI_Op_create(dd_amn_fn, 1, &MPI_DD_AMN);
    MPI_Op_create(zz_amx_fn, 1, &MPI_ZZ_AMX);
    MPI_Op_create(zz_amn_fn, 1, &MPI_ZZ_AMN);

    initialized = 1;
}
