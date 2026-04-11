/* multifloats_c.h -- plain-C companion to the Fortran multifloats module.
 *
 * Provides:
 *   - C-ABI struct types matching the layout of Fortran
 *         type, sequence :: float64x2   real(dp) :: limbs(2)    end type
 *         type, sequence :: complex128x2 type(float64x2) :: re, im  end type
 *   - Inline double-double primitives (Knuth two-sum) used by the migrated
 *     BLACS reduction kernels.
 *   - Global MPI datatype and user-defined MPI_Op handles used by the
 *     migrated ddgsum2d_/ddgamx2d_/ddgamn2d_ routines.
 *   - multifloats_mpi_init(), an idempotent initializer that registers the
 *     MPI handles.  Called from the migrated blacs_gridinit_/blacs_gridmap_
 *     entry points so user code does not need to do anything special.
 */
#ifndef MULTIFLOATS_C_H
#define MULTIFLOATS_C_H 1

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------- */
/* Types                                                                  */
/* --------------------------------------------------------------------- */

typedef struct { double limbs[2]; } float64x2_t;
typedef struct { float64x2_t re, im; } complex128x2_t;

/* --------------------------------------------------------------------- */
/* Scalar primitives                                                      */
/* --------------------------------------------------------------------- */

/* Knuth two-sum: s + e = a + b exactly. */
static inline void mf_two_sum_d(double a, double b, double *s, double *e) {
    double ss = a + b;
    double bb = ss - a;
    *s = ss;
    *e = (a - (ss - bb)) + (b - bb);
}

/* Double-double addition (Briggs / QD add_dd_dd).                        */
static inline float64x2_t mf_add(float64x2_t x, float64x2_t y) {
    double s1, e1, s2, e2;
    mf_two_sum_d(x.limbs[0], y.limbs[0], &s1, &e1);
    mf_two_sum_d(x.limbs[1], y.limbs[1], &s2, &e2);
    double lo = e1 + s2;
    double hi1, lo1;
    mf_two_sum_d(s1, lo, &hi1, &lo1);
    lo1 += e2;
    double hi, loo;
    mf_two_sum_d(hi1, lo1, &hi, &loo);
    float64x2_t r;
    r.limbs[0] = hi;
    r.limbs[1] = loo;
    return r;
}

/* Unary negation. */
static inline float64x2_t mf_neg(float64x2_t x) {
    float64x2_t r;
    r.limbs[0] = -x.limbs[0];
    r.limbs[1] = -x.limbs[1];
    return r;
}

/* Absolute value -- sign of a double-double is determined by the high
 * limb, and negation is distributive so we flip both limbs together. */
static inline float64x2_t mf_abs(float64x2_t x) {
    if (x.limbs[0] < 0.0) {
        x.limbs[0] = -x.limbs[0];
        x.limbs[1] = -x.limbs[1];
    }
    return x;
}

/* Lexicographic compare on (hi, lo). Returns -1, 0, +1. */
static inline int mf_cmp(float64x2_t a, float64x2_t b) {
    if (a.limbs[0] < b.limbs[0]) return -1;
    if (a.limbs[0] > b.limbs[0]) return  1;
    if (a.limbs[1] < b.limbs[1]) return -1;
    if (a.limbs[1] > b.limbs[1]) return  1;
    return 0;
}

/* Complex double-double addition. */
static inline complex128x2_t mf_cadd(complex128x2_t a, complex128x2_t b) {
    complex128x2_t r;
    r.re = mf_add(a.re, b.re);
    r.im = mf_add(a.im, b.im);
    return r;
}

/* |z|_1 = |re| + |im| surrogate used by BLACS Cabs1 in the single/double
 * complex amx/amn kernels. We match that surrogate here so that ddgamx2d
 * and zzgamx2d report the same element that dgamx2d/zgamx2d would have
 * reported in stock BLACS. */
static inline float64x2_t mf_cabs1(complex128x2_t z) {
    return mf_add(mf_abs(z.re), mf_abs(z.im));
}

/* PBLAS quick-return shorthand. PBLAS C entry points test scalar
 * arguments with idioms like ``ALPHA[REAL_PART] == ZERO`` to skip
 * trivial work. C ``==`` is not defined on the float64x2_t struct,
 * so the migrator rewrites those expressions to use these macros.
 * A normalized double-double is zero iff both limbs are bitwise zero;
 * one iff the high limb is exactly 1.0 and the low limb is zero. */
#define MF_IS_ZERO(x) ((x).limbs[0] == 0.0 && (x).limbs[1] == 0.0)
#define MF_IS_ONE(x)  ((x).limbs[0] == 1.0 && (x).limbs[1] == 0.0)

/* --------------------------------------------------------------------- */
/* MPI handles                                                            */
/* --------------------------------------------------------------------- */

extern MPI_Datatype MPI_FLOAT64X2;
extern MPI_Datatype MPI_COMPLEX128X2;

extern MPI_Op MPI_DD_SUM;
extern MPI_Op MPI_ZZ_SUM;
extern MPI_Op MPI_DD_AMX;
extern MPI_Op MPI_DD_AMN;
extern MPI_Op MPI_ZZ_AMX;
extern MPI_Op MPI_ZZ_AMN;

/* Idempotent; safe to call multiple times. Must be called after MPI_Init
 * because it invokes MPI_Type_contiguous / MPI_Op_create. BLACS callers
 * always call blacs_gridinit or blacs_gridmap after MPI_Init, so the
 * migrated entry points trigger this on first use. */
void multifloats_mpi_init(void);

#ifdef __cplusplus
}
#endif

#endif /* MULTIFLOATS_C_H */
