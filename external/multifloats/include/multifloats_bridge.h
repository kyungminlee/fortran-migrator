/* multifloats_bridge.h -- C++ bridge for migrated BLACS / PBLAS sources.
 *
 * Maps the migrated C type names (float64x2_t, cmplxDD, complex128x2_t)
 * to multifloats::MultiFloat<double,2> which provides full operator
 * overloading: ==, !=, <, >, +, -, *, /, +=, -=, unary -, etc.
 *
 * Also declares the runtime MPI handle registration used by the migrated
 * BLACS reductions (MPI_FLOAT64X2, MPI_DD_SUM, etc.) and the
 * multifloats_mpi_init() function that registers them.
 *
 * Usage: compile with  -include multifloats_bridge.h
 *
 * Requires: multifloats.hh (the C++ template library) on the include path.
 */
#ifndef MULTIFLOATS_BRIDGE_H
#define MULTIFLOATS_BRIDGE_H

#ifdef __cplusplus

#include "multifloats.hh"
#include <mpi.h>

/* ------------------------------------------------------------------ */
/* Type aliases matching the names emitted by the C migrator.         */
/* ------------------------------------------------------------------ */

using float64x2_t = multifloats::float64x2;
typedef float64x2_t cmplxDD[2];
struct complex128x2_t { float64x2_t re, im; };

/* ------------------------------------------------------------------ */
/* MF_IS_* macros injected by c_migrator for C scalar comparisons.    */
/* In C++ they just use the native == operator on MultiFloat.         */
/* ------------------------------------------------------------------ */

#define MF_IS_ZERO(x) ((x) == float64x2_t{0.0})
#define MF_IS_ONE(x)  ((x) == float64x2_t{1.0})

/* ------------------------------------------------------------------ */
/* Mixed-type operators.                                              */
/*                                                                    */
/* C++ does not check the right operand for implicit conversion when  */
/* the left operand's type has no matching operator. PBLAS uses       */
/* idioms like ``ONE + ssq * ...`` where ONE=1.0 (double).           */
/* ------------------------------------------------------------------ */

inline float64x2_t operator+(double a, float64x2_t b) { return float64x2_t{a} + b; }
inline float64x2_t operator-(double a, float64x2_t b) { return float64x2_t{a} - b; }
inline float64x2_t operator*(double a, float64x2_t b) { return float64x2_t{a} * b; }
inline float64x2_t operator/(double a, float64x2_t b) { return float64x2_t{a} / b; }
inline bool operator==(double a, float64x2_t b) { return float64x2_t{a} == b; }
inline bool operator!=(double a, float64x2_t b) { return float64x2_t{a} != b; }
inline bool operator<(double a, float64x2_t b)  { return float64x2_t{a} < b; }
inline bool operator>(double a, float64x2_t b)  { return float64x2_t{a} > b; }
inline bool operator<=(double a, float64x2_t b) { return float64x2_t{a} <= b; }
inline bool operator>=(double a, float64x2_t b) { return float64x2_t{a} >= b; }

/* ------------------------------------------------------------------ */
/* Integer truncation helper.                                          */
/* ------------------------------------------------------------------ */

inline int mf_to_int(float64x2_t x) {
    return static_cast<int>(static_cast<double>(x));
}

/* ------------------------------------------------------------------ */
/* Absolute value (for BLACS ABS macro compatibility).                */
/* ------------------------------------------------------------------ */

inline float64x2_t mf_abs(float64x2_t x) { return x < float64x2_t{0.0} ? -x : x; }
inline float64x2_t mf_cabs1(complex128x2_t z) { return mf_abs(z.re) + mf_abs(z.im); }

/* ------------------------------------------------------------------ */
/* MPI handle externs + init (defined in multifloats_mpi.cpp).        */
/* ------------------------------------------------------------------ */

extern "C" {

extern MPI_Datatype MPI_FLOAT64X2;
extern MPI_Datatype MPI_COMPLEX128X2;

extern MPI_Op MPI_DD_SUM;
extern MPI_Op MPI_ZZ_SUM;
extern MPI_Op MPI_DD_AMX;
extern MPI_Op MPI_DD_AMN;
extern MPI_Op MPI_ZZ_AMX;
extern MPI_Op MPI_ZZ_AMN;

void multifloats_mpi_init(void);

} /* extern "C" */

#endif /* __cplusplus */
#endif /* MULTIFLOATS_BRIDGE_H */
