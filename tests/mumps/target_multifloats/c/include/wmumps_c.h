/*
 * Multifloats wrapper around upstream `zmumps_c.h` for the migrated
 * wmumps archive. Mirror of `mmumps_c.h` for the complex side.
 */

#ifndef WMUMPS_C_H
#define WMUMPS_C_H

#ifdef ZMUMPS_C_H
#  error "wmumps_c.h must be included before zmumps_c.h"
#endif

#define ZMUMPS_STRUC_C WMUMPS_STRUC_C
#define zmumps_c       wmumps_c
#define zmumps_f77_    wmumps_f77_

#include "mumps_c_types.h"
#include "zmumps_c.h"

#endif /* WMUMPS_C_H */
