/*
 * Multifloats (double-double / float64x2) override of MUMPS's standard
 * mumps_c_types.h. Mirror of the kind16 / kind10 shadows, using a 16-byte
 * POD struct to model upstream multifloats's ``float64x2`` (two
 * double-precision limbs) for the real width, and 32-byte struct for
 * complex (real-imaginary pair of float64x2).
 *
 * Layout matters: the migrated mmumps / wmumps Fortran archives use
 * ``TYPE(real64x2)`` / ``TYPE(cmplx64x2)`` which gfortran lays out as a
 * sequence of two / four ``REAL(8)`` words. Our C-side struct mirrors
 * that — both sides see 16 / 32 contiguous bytes. The C bridge passes
 * these by pointer through to the Fortran entry; no value-level
 * arithmetic is done on the C side, so a plain POD shadow is enough.
 *
 * The full operator-rich ``multifloats::float64x2`` C++ class lives in
 * external/multifloats-mpi/multifloats_bridge.h; that header is C++
 * only and not appropriate for mumps_c.c (which is plain C).
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

#define SMUMPS_COMPLEX float
#define SMUMPS_REAL    float

/* Plain-C POD layouts for multifloats real / complex. The migrator's
 * Fortran TYPE(real64x2) is two REAL(8); TYPE(cmplx64x2) is the
 * obvious .re/.im pair, total four REAL(8). Match byte-for-byte. */
typedef struct { double limbs[2]; }      mumps_float64x2;
typedef struct { mumps_float64x2 r, i; } mumps_complex64x2;

#define DMUMPS_COMPLEX mumps_float64x2
#define DMUMPS_REAL    mumps_float64x2

typedef struct { float       r, i; } mumps_complex;
typedef mumps_complex64x2            mumps_double_complex;

#define CMUMPS_COMPLEX mumps_complex
#define CMUMPS_REAL    float

#define ZMUMPS_COMPLEX mumps_complex64x2
#define ZMUMPS_REAL    mumps_float64x2


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
