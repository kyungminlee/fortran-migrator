/*
 * Basic smoke test for the XMUMPS C bridge — mirror of
 * test_dmumps_c_basic.c with complex (__float128 r/i) entries.
 *
 * mumps_double_complex is `mumps_double_complex` after our shadow
 * mumps_c_types.h overrides it to `{ __float128 r, i; }` — gfortran's
 * COMPLEX(KIND=16) ABI is two consecutive __float128 values, so a
 * struct of two __float128 fields with no extra padding is
 * binary-compatible.
 */

#include <math.h>
#include <quadmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include "xmumps_c.h"

static FILE *gJson = NULL;
static int gAnyFail = 0, gCaseCount = 0;

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
    printf("  test_zmumps_c_basic [%s] max_rel_err=%s  %s\n",
           case_label, relbuf, passed ? "PASS" : "FAIL");
}

static void report_finalize_c(void)
{
    fprintf(gJson, "  ]\n}\n");
    fclose(gJson);
}

static int report_status_c(void) { return gAnyFail; }

/* Multiplication of mumps_double_complex = {r, i}. */
static mumps_double_complex cmul(mumps_double_complex a, mumps_double_complex b)
{
    mumps_double_complex c;
    c.r = a.r * b.r - a.i * b.i;
    c.i = a.r * b.i + a.i * b.r;
    return c;
}

static __float128 cabs_q(mumps_double_complex z)
{
    return sqrtq(z.r * z.r + z.i * z.i);
}

int main(int argc, char **argv)
{
    XMUMPS_STRUC_C id = {0};   /* L-5 — see test_dmumps_c_basic.c */
    enum { N = 4 };

    /* Diagonally-dominant complex unsymmetric A. Hand-built so the
     * test is reproducible and the residual stays bounded. */
    mumps_double_complex A[N][N] = {
        {{ 5.0q,  0.5q}, { 0.5q,  0.1q}, {-0.25q, 0.05q}, { 0.1q,  0.0q}},
        {{ 0.3q,  0.0q}, {-6.0q, -0.4q}, { 0.5q,  0.1q},  { 0.2q, -0.1q}},
        {{-0.4q,  0.1q}, { 0.2q, -0.05q},{ 7.0q,  0.3q},  {-0.3q,  0.2q}},
        {{ 0.1q, -0.05q},{-0.3q,  0.1q}, { 0.4q,  0.0q},  {-8.0q,  0.5q}},
    };
    mumps_double_complex x_true[N] = {
        { 1.0q, -0.5q}, {-2.0q,  1.0q}, { 3.0q, -1.5q}, {-4.0q,  2.0q}
    };

    MUMPS_INT       irn[N*N], jcn[N*N];
    mumps_double_complex a_vals[N*N], rhs[N];

    MPI_Init(&argc, &argv);
    {
        int i, j, k = 0;
        for (j = 0; j < N; j++)
            for (i = 0; i < N; i++, k++) {
                irn[k] = i+1; jcn[k] = j+1; a_vals[k] = A[i][j];
            }
        /* b = A * x_true */
        for (i = 0; i < N; i++) {
            mumps_double_complex s = { 0.0q, 0.0q };
            for (j = 0; j < N; j++) {
                mumps_double_complex p = cmul(A[i][j], x_true[j]);
                s.r += p.r; s.i += p.i;
            }
            rhs[i] = s;
        }
    }

    report_init_c("test_zmumps_c_basic", "kind16");

    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    xmumps_c(&id);
    if (id.infog[0] < 0) { fprintf(stderr, "JOB=-1 failed: %d\n", id.infog[0]); return 1; }

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    id.n   = N;
    id.nnz = (MUMPS_INT8) (N * N);
    id.irn = irn;
    id.jcn = jcn;
    id.a   = a_vals;
    id.rhs = rhs;

    id.job = 6;
    xmumps_c(&id);
    if (id.infog[0] < 0) {
        fprintf(stderr, "JOB=6 failed: infog=%d %d\n", id.infog[0], id.infog[1]);
        return 1;
    }

    {
        __float128 max_rel = 0.0q, denom = 0.0q;
        int i;
        for (i = 0; i < N; i++) {
            __float128 a = cabs_q(x_true[i]);
            if (a > denom) denom = a;
        }
        for (i = 0; i < N; i++) {
            mumps_double_complex d = { rhs[i].r - x_true[i].r, rhs[i].i - x_true[i].i };
            __float128 m = cabs_q(d);
            if (m > max_rel) max_rel = m;
        }
        if (denom > FLT128_MIN) max_rel /= denom;
        report_case_c("n=4", max_rel,
                      16.0q * (__float128)(N*N*N) * FLT128_EPSILON);
    }

    id.job = -2;
    xmumps_c(&id);

    report_finalize_c();
    MPI_Finalize();
    return report_status_c() ? 1 : 0;
}
