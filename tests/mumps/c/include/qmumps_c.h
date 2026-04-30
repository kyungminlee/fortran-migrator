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

/* CRITICAL: include OUR mumps_c_types.h shadow FIRST so its
 * `#define MUMPS_C_TYPES_H` guard fires. dmumps_c.h's own
 * `#include "mumps_c_types.h"` then resolves to a no-op (guard
 * already set), which means dmumps_c.h instantiates DMUMPS_STRUC_C
 * with our __float128 typedefs instead of upstream's `double`.
 *
 * Without this prelude, dmumps_c.h's `#include "mumps_c_types.h"`
 * would search dmumps_c.h's own directory first (== upstream's
 * include/) and find upstream's `double`-typed mumps_c_types.h
 * before our shadow on the -I path. */
#include "mumps_c_types.h"
#include "dmumps_c.h"

#endif /* QMUMPS_C_H */
