/*
 * Basic smoke test for the QMUMPS C bridge (dmumps_c → qmumps_c via
 * the include/qmumps_c.h header overrides). Mirrors
 * tests/mumps/fortran/test_dmumps_basic.f90 in shape but operates
 * end-to-end through the C interface.
 *
 * Hand-coded (no quad-precision random helpers shared from Fortran)
 * with a small fixed problem: n=4 diagonally-dominant unsymmetric
 * matrix, x_true known, b = A·x_true. Pass/fail decided by JSON
 * report comparing |x_solve − x_true| / |x_true| against an
 * O(n^3) tolerance scaled by target eps (= __float128 epsilon).
 */

#include <math.h>
#include <quadmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include "qmumps_c.h"

/* Same JSON report layout as prec_report.f90 produces, hand-rolled
 * here so the C test side doesn't need to call back into Fortran. */
static FILE *gJson = NULL;
static int gAnyFail = 0;
static int gCaseCount = 0;

static __float128 target_eps_q(void) { return FLT128_EPSILON; }

static void report_init_c(const char *test_name, const char *target_name)
{
    char filename[256];
    snprintf(filename, sizeof filename, "%s.%s.json", test_name, target_name);
    gJson = fopen(filename, "w");
    if (!gJson) { fprintf(stderr, "report_init_c: cannot open %s\n", filename); exit(1); }
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
            "    {\n"
            "      \"case\":        \"%s\",\n"
            "      \"max_rel_err\": %s,\n"
            "      \"tolerance\":   %s,\n"
            "      \"passed\":      %s\n"
            "    }\n",
            case_label, relbuf, tolbuf, passed ? "true" : "false");
    printf("  %s [%s] max_rel_err=%s  %s\n", "test_dmumps_c_basic",
           case_label, relbuf, passed ? "PASS" : "FAIL");
}

static void report_finalize_c(void)
{
    fprintf(gJson, "  ]\n}\n");
    fclose(gJson);
    gJson = NULL;
    if (gAnyFail) exit(1);
}

int main(int argc, char **argv)
{
    QMUMPS_STRUC_C id;
    int n = 4;
    /* Diagonally-dominant unsymmetric A (column-major triplet). */
    MUMPS_INT  irn[16];
    MUMPS_INT  jcn[16];
    __float128 a_vals[16];
    __float128 x_true[4] = {  1.0q, -2.0q,  3.0q,  -4.0q };
    __float128 rhs[4];

    MPI_Init(&argc, &argv);

    /* Fill A (4x4): diagonal entries large, off-diagonals small.
     * Layout matches dense_to_triplet (column-major flatten). */
    {
        __float128 A[4][4] = {
            {  5.0q,  0.5q, -0.25q,  0.1q  },
            {  0.3q, -6.0q,  0.5q,   0.2q  },
            { -0.4q,  0.2q,  7.0q,  -0.3q  },
            {  0.1q, -0.3q,  0.4q,  -8.0q  },
        };
        int i, j, k = 0;
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++, k++) {
                irn[k]    = i + 1;
                jcn[k]    = j + 1;
                a_vals[k] = A[i][j];
            }
        }
        /* b = A * x_true */
        for (i = 0; i < n; i++) {
            __float128 s = 0.0q;
            for (j = 0; j < n; j++) s += A[i][j] * x_true[j];
            rhs[i] = s;
        }
    }

    report_init_c("test_dmumps_c_basic", "kind16");

    /* ── Init ─────────────────────────────────────────────────────── */
    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    qmumps_c(&id);
    if (id.infog[0] < 0) {
        fprintf(stderr, "JOB=-1 failed: infog[0]=%d\n", id.infog[0]);
        return 1;
    }

    /* Silence diagnostic output. */
    id.icntl[0] = -1;  /* error stream */
    id.icntl[1] = -1;  /* diagnostic */
    id.icntl[2] = -1;  /* global info */
    id.icntl[3] =  0;  /* print level */

    /* ── Problem data ────────────────────────────────────────────── */
    id.n    = n;
    id.nnz  = (MUMPS_INT8) (n * n);
    id.irn  = irn;
    id.jcn  = jcn;
    id.a    = a_vals;
    id.rhs  = rhs;

    /* ── Solve ───────────────────────────────────────────────────── */
    id.job = 6;
    qmumps_c(&id);
    if (id.infog[0] < 0) {
        fprintf(stderr, "JOB=6 failed: infog[0]=%d, infog[1]=%d\n",
                id.infog[0], id.infog[1]);
        return 1;
    }

    /* On exit, rhs[] holds the solution. */
    {
        __float128 max_rel = 0.0q, denom = 0.0q;
        int i;
        for (i = 0; i < n; i++) {
            __float128 ax = fabsq(x_true[i]);
            if (ax > denom) denom = ax;
        }
        for (i = 0; i < n; i++) {
            __float128 d = fabsq(rhs[i] - x_true[i]);
            if (d > max_rel) max_rel = d;
        }
        if (denom > FLT128_MIN) max_rel /= denom;

        __float128 tol = 16.0q * (__float128) (n * n * n) * target_eps_q();
        report_case_c("n=4", max_rel, tol);
    }

    /* ── Cleanup ─────────────────────────────────────────────────── */
    id.job = -2;
    qmumps_c(&id);

    report_finalize_c();
    MPI_Finalize();
    return 0;
}
