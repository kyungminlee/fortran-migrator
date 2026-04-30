/*
 * SYM coverage via the QMUMPS C bridge — mirrors
 * tests/mumps/fortran/test_dmumps_sym.f90 with hand-coded n=4
 * problems for each SYM ∈ {0, 1, 2}.
 *
 * SYM=0  :  general unsymmetric A (LU)
 * SYM=1  :  symmetric positive-definite A = X*X^T + n*I (Cholesky)
 * SYM=2  :  general symmetric A = (R + R^T)/2 + n*I (Bunch-Kaufman)
 *
 * For SYM>0 only the upper triangle is sent in the triplet form; the
 * solver infers the lower triangle.
 */

#include <math.h>
#include <quadmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include "qmumps_c.h"

static FILE *gJson = NULL;
static int gAnyFail = 0;
static int gCaseCount = 0;

static void report_init_c(const char *test_name, const char *target_name)
{
    char filename[256];
    snprintf(filename, sizeof filename, "%s.%s.json", test_name, target_name);
    gJson = fopen(filename, "w");
    if (!gJson) { fprintf(stderr, "report_init: cannot open %s\n", filename); exit(1); }
    fprintf(gJson, "{\n  \"routine\": \"%s\",\n  \"target\":  \"%s\",\n  \"cases\": [\n",
            test_name, target_name);
}

static void report_case_c(const char *case_label, __float128 max_rel, __float128 tol)
{
    char relbuf[64], tolbuf[64];
    int passed = (max_rel <= tol);
    if (!passed) gAnyFail = 1;
    quadmath_snprintf(relbuf, sizeof relbuf, "%.6Qe", max_rel);
    quadmath_snprintf(tolbuf, sizeof tolbuf, "%.6Qe", tol);
    if (gCaseCount > 0) fprintf(gJson, "    ,\n");
    gCaseCount++;
    fprintf(gJson,
            "    {\n      \"case\":        \"%s\",\n"
            "      \"max_rel_err\": %s,\n      \"tolerance\":   %s,\n"
            "      \"passed\":      %s\n    }\n",
            case_label, relbuf, tolbuf, passed ? "true" : "false");
    printf("  test_dmumps_c_sym [%s] max_rel_err=%s  %s\n",
           case_label, relbuf, passed ? "PASS" : "FAIL");
}

static void report_finalize_c(void)
{
    fprintf(gJson, "  ]\n}\n");
    fclose(gJson);
}

static int report_status_c(void) { return gAnyFail; }

/* Solve A*x = b via QMUMPS, given a triplet (irn, jcn, a) of nz entries
 * and an n-vector b (overwritten with x on exit). */
static int qmumps_solve(int sym, MUMPS_INT n, MUMPS_INT8 nnz,
                         MUMPS_INT *irn, MUMPS_INT *jcn,
                         __float128 *a_vals, __float128 *rhs)
{
    QMUMPS_STRUC_C id;
    id.par = 1;
    id.sym = sym;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    qmumps_c(&id);
    if (id.infog[0] < 0) return id.infog[0];

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    id.n   = n;
    id.nnz = nnz;
    id.irn = irn;
    id.jcn = jcn;
    id.a   = a_vals;
    id.rhs = rhs;
    id.job = 6;
    qmumps_c(&id);
    int code = id.infog[0];

    id.job = -2;
    qmumps_c(&id);
    return code;
}

int main(int argc, char **argv)
{
    enum { N = 4 };
    __float128 x_true[N] = { 1.0q, -2.0q, 3.0q, -4.0q };

    MPI_Init(&argc, &argv);
    report_init_c("test_dmumps_c_sym", "kind16");

    /* ── SYM=0 — general unsymmetric ─────────────────────────────── */
    {
        __float128 A[N][N] = {
            {  5.0q,  0.5q, -0.25q,  0.1q  },
            {  0.3q, -6.0q,  0.5q,   0.2q  },
            { -0.4q,  0.2q,  7.0q,  -0.3q  },
            {  0.1q, -0.3q,  0.4q,  -8.0q  },
        };
        MUMPS_INT  irn[N*N], jcn[N*N];
        __float128 vals[N*N], rhs[N];
        int i, j, k = 0;
        for (j = 0; j < N; j++)
            for (i = 0; i < N; i++, k++) {
                irn[k] = i+1; jcn[k] = j+1; vals[k] = A[i][j];
            }
        for (i = 0; i < N; i++) {
            __float128 s = 0.0q;
            for (j = 0; j < N; j++) s += A[i][j] * x_true[j];
            rhs[i] = s;
        }
        int code = qmumps_solve(0, N, (MUMPS_INT8)(N*N), irn, jcn, vals, rhs);
        if (code < 0) { fprintf(stderr, "SYM=0 failed: %d\n", code); return 1; }
        __float128 max_rel = 0.0q, denom = 0.0q;
        for (i = 0; i < N; i++) { __float128 a = fabsq(x_true[i]); if (a > denom) denom = a; }
        for (i = 0; i < N; i++) { __float128 d = fabsq(rhs[i]-x_true[i]); if (d > max_rel) max_rel = d; }
        if (denom > FLT128_MIN) max_rel /= denom;
        report_case_c("sym=0", max_rel,
                      16.0q * (__float128)(N*N*N) * FLT128_EPSILON);
    }

    /* ── SYM=1 — symmetric positive-definite ─────────────────────── */
    {
        /* A = X*X^T + N*I — guaranteed SPD. Hand-build it. */
        __float128 X[N][N] = {
            {  1.0q,  0.5q, -0.25q,  0.1q  },
            {  0.3q,  0.8q,  0.5q,   0.2q  },
            { -0.4q,  0.2q,  0.7q,  -0.3q  },
            {  0.1q, -0.3q,  0.4q,   0.6q  },
        };
        __float128 A[N][N] = { 0 };
        int i, j, k;
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++) {
                __float128 s = 0.0q;
                for (k = 0; k < N; k++) s += X[i][k] * X[j][k];
                A[i][j] = s + (i == j ? (__float128) N : 0.0q);
            }
        /* Upper-triangle triplet for SYM=1 */
        MUMPS_INT  irn[N*(N+1)/2], jcn[N*(N+1)/2];
        __float128 vals[N*(N+1)/2], rhs[N];
        int idx = 0;
        for (j = 0; j < N; j++)
            for (i = 0; i <= j; i++, idx++) {
                irn[idx] = i+1; jcn[idx] = j+1; vals[idx] = A[i][j];
            }
        for (i = 0; i < N; i++) {
            __float128 s = 0.0q;
            for (j = 0; j < N; j++) s += A[i][j] * x_true[j];
            rhs[i] = s;
        }
        int code = qmumps_solve(1, N, (MUMPS_INT8) idx, irn, jcn, vals, rhs);
        if (code < 0) { fprintf(stderr, "SYM=1 failed: %d\n", code); return 1; }
        __float128 max_rel = 0.0q, denom = 0.0q;
        for (i = 0; i < N; i++) { __float128 a = fabsq(x_true[i]); if (a > denom) denom = a; }
        for (i = 0; i < N; i++) { __float128 d = fabsq(rhs[i]-x_true[i]); if (d > max_rel) max_rel = d; }
        if (denom > FLT128_MIN) max_rel /= denom;
        report_case_c("sym=1", max_rel,
                      16.0q * (__float128)(N*N*N) * FLT128_EPSILON);
    }

    /* ── SYM=2 — general symmetric ───────────────────────────────── */
    {
        __float128 A[N][N] = {
            {  6.0q,  0.5q, -0.25q,  0.1q  },
            {  0.5q, -7.0q,  0.5q,   0.2q  },
            { -0.25q, 0.5q,  8.0q,  -0.3q  },
            {  0.1q,  0.2q, -0.3q,  -9.0q  },
        };  /* hand-symmetric */
        MUMPS_INT  irn[N*(N+1)/2], jcn[N*(N+1)/2];
        __float128 vals[N*(N+1)/2], rhs[N];
        int i, j, idx = 0;
        for (j = 0; j < N; j++)
            for (i = 0; i <= j; i++, idx++) {
                irn[idx] = i+1; jcn[idx] = j+1; vals[idx] = A[i][j];
            }
        for (i = 0; i < N; i++) {
            __float128 s = 0.0q;
            for (j = 0; j < N; j++) s += A[i][j] * x_true[j];
            rhs[i] = s;
        }
        int code = qmumps_solve(2, N, (MUMPS_INT8) idx, irn, jcn, vals, rhs);
        if (code < 0) { fprintf(stderr, "SYM=2 failed: %d\n", code); return 1; }
        __float128 max_rel = 0.0q, denom = 0.0q;
        for (i = 0; i < N; i++) { __float128 a = fabsq(x_true[i]); if (a > denom) denom = a; }
        for (i = 0; i < N; i++) { __float128 d = fabsq(rhs[i]-x_true[i]); if (d > max_rel) max_rel = d; }
        if (denom > FLT128_MIN) max_rel /= denom;
        report_case_c("sym=2", max_rel,
                      16.0q * (__float128)(N*N*N) * FLT128_EPSILON);
    }

    report_finalize_c();
    MPI_Finalize();
    return report_status_c() ? 1 : 0;
}
