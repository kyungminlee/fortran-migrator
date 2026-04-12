/* pblas_mf_bridge.h -- C++ bridge between migrated PBLAS sources and
 * the multifloats C++ template library (multifloats.hh).
 *
 * When PBLAS .c files are compiled as C++ (mpicxx -std=c++17), this
 * header maps the migrated C type names (float64x2_t, cmplxDD,
 * complex128x2_t) to multifloats::MultiFloat<double,2> which provides
 * full operator overloading: ==, !=, <, >, +, -, *, /, +=, -=, etc.
 * This eliminates most struct-operator issues that plague the pure-C
 * compilation path.
 *
 * Usage: compile with  -include pblas_mf_bridge.h  -DMULTIFLOATS_C_H=1
 * (the guard macro suppresses the plain-C multifloats_c.h which would
 * conflict with the C++ type aliases defined here).
 */
#ifndef PBLAS_MF_BRIDGE_H
#define PBLAS_MF_BRIDGE_H

#ifdef __cplusplus
#include "multifloats.hh"

using float64x2_t = multifloats::float64x2;
typedef float64x2_t cmplxDD[2];
struct complex128x2_t { float64x2_t re, im; };

/* ------------------------------------------------------------------ */
/* MF_IS_* macros injected by c_migrator for C scalar comparisons.    */
/* In C++ they just use the native == operator on MultiFloat.         */
/* ------------------------------------------------------------------ */
#define MF_IS_ZERO(x)  ((x) == float64x2_t{0.0})
#define MF_IS_ONE(x)   ((x) == float64x2_t{1.0})

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
/* Integer truncation.                                                */
/*                                                                    */
/* PBLAS packs integer indices into work buffers as double scalars;   */
/* (Int)(work[1]) in the original source becomes                      */
/* (Int)(float64x2_t) in the migrated code. C++ explicit operator T() */
/* blocks the direct cast, so the c_migrator rewrites the pattern     */
/* into mf_to_int(work[1]).                                           */
/* ------------------------------------------------------------------ */
inline int mf_to_int(float64x2_t x) {
    return static_cast<int>(static_cast<double>(x));
}

#endif /* __cplusplus */
#endif /* PBLAS_MF_BRIDGE_H */
