/*
 * Quad-precision (kind16) wrapper around upstream `dmumps_c.h`.
 *
 * Use this instead of `dmumps_c.h` from C tests / drivers that target the
 * migrated qmumps archive. It macro-renames the entry point (`dmumps_c`
 * → `qmumps_c`), the F77 symbol (`dmumps_f77_` → `qmumps_f77_`) and the
 * struct typename (`DMUMPS_STRUC_C` → `QMUMPS_STRUC_C`), then includes
 * the upstream header — which itself pulls in our quad-precision
 * `mumps_c_types.h` shadow and instantiates the struct with
 * `__float128` numeric fields.
 *
 * The same .o file produced by compiling upstream `mumps_c.c` with
 *   -DMUMPS_ARITH=MUMPS_ARITH_d -Ddmumps_f77_=qmumps_f77_
 *   -Ddmumps_c=qmumps_c -DDMUMPS_STRUC_C=QMUMPS_STRUC_C
 * provides the symbols `qmumps_c` (entry) and the static state for
 * `dmumps_assign_*` / `dmumps_nullify_c_*` (kept under the original
 * D prefix because that is what the migrated qmumps_f77_ expects to
 * call back into).
 */

#ifndef QMUMPS_C_H
#define QMUMPS_C_H

#define DMUMPS_STRUC_C QMUMPS_STRUC_C
#define dmumps_c       qmumps_c
#define dmumps_f77_    qmumps_f77_

#include "dmumps_c.h"

#endif /* QMUMPS_C_H */
