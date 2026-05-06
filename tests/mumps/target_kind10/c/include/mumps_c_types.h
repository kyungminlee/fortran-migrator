/*
 * Extended-precision (kind10 = x87 ``long double``) override of MUMPS's
 * standard mumps_c_types.h. Mirror of the kind16 shadow at
 * tests/mumps/target_kind16/c/include/mumps_c_types.h, switched to
 * 80-bit ``long double`` for the real / complex widths the migrated
 * emumps / ymumps archives expect (REAL(KIND=10) / COMPLEX(KIND=10)
 * in gfortran).
 *
 * Selected via include-path priority — when this file's directory is
 * on the include path BEFORE the upstream MUMPS ``include/``, every
 * consumer of ``<mumps_c_types.h>`` sees the kind10 typedefs instead
 * of ``double``.
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

/* The two REAL/COMPLEX widths the migrated emumps / ymumps archives
 * actually expect — REAL(KIND=10) and COMPLEX(KIND=10) in gfortran
 * (== ``long double`` on x86-64 with the Intel ABI 80-bit format,
 * sizeof 16 bytes after alignment). */
#define DMUMPS_COMPLEX long double
#define DMUMPS_REAL    long double

typedef struct { float       r, i; } mumps_complex;
typedef struct { long double r, i; } mumps_double_complex;

#define CMUMPS_COMPLEX mumps_complex
#define CMUMPS_REAL    float

#define ZMUMPS_COMPLEX mumps_double_complex
#define ZMUMPS_REAL    long double


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
