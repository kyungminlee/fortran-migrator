/*
 * Multifloats wrapper around upstream `dmumps_c.h` for the migrated
 * mmumps archive. Mirror of target_kind16/c/include/qmumps_c.h with
 * the arithmetic-prefix renamed to m/M.
 */

#ifndef MMUMPS_C_H
#define MMUMPS_C_H

#ifdef DMUMPS_C_H
#  error "mmumps_c.h must be included before dmumps_c.h"
#endif

#define DMUMPS_STRUC_C MMUMPS_STRUC_C
#define dmumps_c       mmumps_c
#define dmumps_f77_    mmumps_f77_

#include "mumps_c_types.h"
#include "dmumps_c.h"

#endif /* MMUMPS_C_H */
