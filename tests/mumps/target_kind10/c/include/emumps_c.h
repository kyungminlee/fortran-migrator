/*
 * kind10 wrapper around upstream `dmumps_c.h` for the migrated emumps
 * archive. Mirror of target_kind16/c/include/qmumps_c.h with the
 * arithmetic-prefix renamed to e/E.
 *
 * Use this instead of `dmumps_c.h` from C tests / drivers that target
 * the migrated emumps archive. It macro-renames the entry point
 * (`dmumps_c` → `emumps_c`), the F77 symbol (`dmumps_f77_` →
 * `emumps_f77_`), and the struct typename (`DMUMPS_STRUC_C` →
 * `EMUMPS_STRUC_C`), then includes the upstream header — which itself
 * pulls in our kind10 `mumps_c_types.h` shadow and instantiates the
 * struct with `long double` numeric fields.
 */

#ifndef EMUMPS_C_H
#define EMUMPS_C_H

/* L-1 guard: detect wrong include order. Including upstream's
 * dmumps_c.h FIRST sets DMUMPS_C_H and instantiates DMUMPS_STRUC_C
 * BEFORE our renames have a chance to redirect the typename — a
 * later include of this header would silently use the un-renamed
 * struct. Fail loudly instead. */
#ifdef DMUMPS_C_H
#  error "emumps_c.h must be included before dmumps_c.h"
#endif

#define DMUMPS_STRUC_C EMUMPS_STRUC_C
#define dmumps_c       emumps_c
#define dmumps_f77_    emumps_f77_

#include "mumps_c_types.h"
#include "dmumps_c.h"

#endif /* EMUMPS_C_H */
