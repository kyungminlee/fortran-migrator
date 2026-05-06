/*
 * kind10 wrapper around upstream `zmumps_c.h` for the migrated ymumps
 * archive. Mirror of `emumps_c.h` for the complex-arithmetic side; see
 * that file for design rationale.
 */

#ifndef YMUMPS_C_H
#define YMUMPS_C_H

#ifdef ZMUMPS_C_H
#  error "ymumps_c.h must be included before zmumps_c.h"
#endif

#define ZMUMPS_STRUC_C YMUMPS_STRUC_C
#define zmumps_c       ymumps_c
#define zmumps_f77_    ymumps_f77_

#include "mumps_c_types.h"
#include "zmumps_c.h"

#endif /* YMUMPS_C_H */
