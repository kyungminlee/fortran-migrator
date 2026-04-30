/*
 * Quad-precision (kind16) override of MUMPS's standard mumps_c_types.h.
 *
 * Shadows external/MUMPS_5.8.2/include/mumps_c_types.h via include-path
 * priority — when this file is on the include path BEFORE the upstream
 * `include/` directory, every consumer of `<mumps_c_types.h>` (upstream
 * `mumps_c.c`, `mumps_common.c`, `mumps_addr.c`, `dmumps_c.h`,
 * `zmumps_c.h`, ...) sees the quad typedefs instead of `double`.
 *
 * Match upstream's symbol surface byte-for-byte except for the
 * REAL/COMPLEX widths, so upstream `mumps_c.c` compiles unchanged.
 */

#ifndef MUMPS_C_TYPES_H
#define MUMPS_C_TYPES_H

#include <stdint.h>

#include "mumps_int_def.h"

#ifdef MUMPS_INTSIZE64
#define MUMPS_INT int64_t
#else
#define MUMPS_INT int
#endif

#define MUMPS_INT8 int64_t

/* SMUMPS / CMUMPS unchanged — bridge does not expose single-precision
 * arithmetics, but keeping the macros lets upstream headers compile. */
#define SMUMPS_COMPLEX float
#define SMUMPS_REAL    float

/* The two REAL/COMPLEX widths the migrated qmumps / xmumps archives
 * actually expect — REAL(KIND=16) and COMPLEX(KIND=16) in gfortran. */
#define DMUMPS_COMPLEX __float128
#define DMUMPS_REAL    __float128

typedef struct { float       r, i; } mumps_complex;
typedef struct { __float128  r, i; } mumps_double_complex;

#define CMUMPS_COMPLEX mumps_complex
#define CMUMPS_REAL    float

#define ZMUMPS_COMPLEX mumps_double_complex
#define ZMUMPS_REAL    __float128


#ifndef mumps_ftnlen
# define mumps_ftnlen MUMPS_INT
#endif


#define MUMPS_ARITH_s 1
#define MUMPS_ARITH_d 2
#define MUMPS_ARITH_c 4
#define MUMPS_ARITH_z 8

#define MUMPS_ARITH_REAL   ( MUMPS_ARITH_s | MUMPS_ARITH_d )
#define MUMPS_ARITH_CMPLX  ( MUMPS_ARITH_c | MUMPS_ARITH_z )
#define MUMPS_ARITH_SINGLE ( MUMPS_ARITH_s | MUMPS_ARITH_c )
#define MUMPS_ARITH_DBL    ( MUMPS_ARITH_d | MUMPS_ARITH_z )

#define MUMPS_OFF_T MUMPS_INT8

#endif /* MUMPS_C_TYPES_H */
