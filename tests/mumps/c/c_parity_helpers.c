/*
 * Fortran-callable helpers that drive the C bridge (qmumps_c /
 * xmumps_c). They exist solely for the parity tests: a Fortran
 * driver can invoke `c_qmumps_solve` / `c_xmumps_solve` to run a
 * solve through the C bridge and compare the result bit-for-bit
 * against the same problem solved via target_qmumps / target_xmumps
 * (Fortran path).
 *
 * Both paths converge to the same migrated qmumps_f77_ / xmumps_f77_
 * Fortran entry; if the C struct extraction in mumps_c.c diverges
 * from the Fortran-direct path, the parity test surfaces it as a
 * non-zero diff.
 */

#include "qmumps_c.h"
#include "xmumps_c.h"
#include <mpi.h>

int c_qmumps_solve(int *n_in, MUMPS_INT8 *nnz_in,
                   MUMPS_INT *irn, MUMPS_INT *jcn,
                   __float128 *a_vals, __float128 *rhs)
{
    QMUMPS_STRUC_C id = {0};
    int code;

    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    qmumps_c(&id);
    if (id.infog[0] < 0) return id.infog[0];

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    id.n   = *n_in;
    id.nnz = *nnz_in;
    id.irn = irn;
    id.jcn = jcn;
    id.a   = a_vals;
    id.rhs = rhs;
    id.job = 6;
    qmumps_c(&id);
    code = id.infog[0];

    id.job = -2;
    qmumps_c(&id);
    return code;
}

int c_xmumps_solve(int *n_in, MUMPS_INT8 *nnz_in,
                   MUMPS_INT *irn, MUMPS_INT *jcn,
                   mumps_double_complex *a_vals, mumps_double_complex *rhs)
{
    XMUMPS_STRUC_C id = {0};
    int code;

    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    xmumps_c(&id);
    if (id.infog[0] < 0) return id.infog[0];

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    id.n   = *n_in;
    id.nnz = *nnz_in;
    id.irn = irn;
    id.jcn = jcn;
    id.a   = a_vals;
    id.rhs = rhs;
    id.job = 6;
    xmumps_c(&id);
    code = id.infog[0];

    id.job = -2;
    xmumps_c(&id);
    return code;
}
