/*
 * Quad-precision (kind16) wrapper around upstream `zmumps_c.h`.
 *
 * Mirror of `qmumps_c.h` for the complex-arithmetic side. See that file
 * for design rationale. The matching object file is built by compiling
 * upstream `mumps_c.c` with
 *   -DMUMPS_ARITH=MUMPS_ARITH_z -Dzmumps_f77_=xmumps_f77_
 *   -Dzmumps_c=xmumps_c -DZMUMPS_STRUC_C=XMUMPS_STRUC_C
 * which exposes the entry point `xmumps_c` plus static state for
 * `zmumps_assign_*` / `zmumps_nullify_c_*` (kept under Z because that
 * is what the migrated xmumps_f77_ calls back into).
 */

#ifndef XMUMPS_C_H
#define XMUMPS_C_H

#define ZMUMPS_STRUC_C XMUMPS_STRUC_C
#define zmumps_c       xmumps_c
#define zmumps_f77_    xmumps_f77_

#include "zmumps_c.h"

#endif /* XMUMPS_C_H */
