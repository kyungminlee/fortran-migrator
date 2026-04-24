/* multifloats_bridge.h -- C++ bridge for migrated BLACS / PBLAS sources.
 *
 * Maps the migrated C type names (float64x2_t, cmplxDD, complex128x2_t)
 * to multifloats::float64x2, the operator-rich C++ class from the
 * upstream library (github.com/kyungminlee/multifloats). Provides:
 *
 *   - full operator overloading (==, !=, <, >, +, -, *, /, +=, -=, ...)
 *   - mixed-type operators (``double op float64x2_t``)
 *   - MF_IS_ZERO / MF_IS_ONE macros for the C migrator's PBLAS rewrites
 *   - MPI datatype + op handle declarations (registered by
 *     multifloats_mpi_init, defined in multifloats_mpi.cpp)
 *
 * Usage: compile with  -include multifloats_bridge.h
 *
 * Requires: the multifloats headers on the include path (either
 * installed, or from a FetchContent'd checkout).
 */
#ifndef MULTIFLOATS_BRIDGE_H
#define MULTIFLOATS_BRIDGE_H

#ifdef __cplusplus

/* Pull in the upstream multifloats header, but hide its plain-C POD
 * ``float64x2_t`` typedef under an internal name so our own
 * ``float64x2_t`` (the operator-rich ``multifloats::float64x2`` alias
 * below) can reuse the public identifier. The two types have identical
 * layout (the upstream header static_asserts it). Migrated C-as-C++
 * code never touches the POD form — the ``extern "C"`` dd_* prototypes
 * that use it are not called from migrated source.
 *
 * PBLAS's ``pblas.h`` defines ``#define Int int`` at file scope, which
 * breaks multifloats's ``template <typename Int>`` parameter name.
 * Push/pop the macro around the upstream include so the template
 * parameter is parsed verbatim while leaving the rest of the TU's use
 * of ``Int`` alone.
 */
#pragma push_macro("Int")
#undef Int
#define float64x2_t   mf_upstream_float64x2_pod_t
#define complex64x2_t mf_upstream_complex64x2_pod_t
#include <multifloats.h>
#undef float64x2_t
#undef complex64x2_t
#pragma pop_macro("Int")

#include <mpi.h>

/* ------------------------------------------------------------------ */
/* Type aliases matching the names emitted by the C migrator.         */
/* ------------------------------------------------------------------ */

using float64x2_t = multifloats::float64x2;
typedef float64x2_t cmplxDD[2];
struct complex128x2_t {
    union { float64x2_t r; float64x2_t re; };
    union { float64x2_t i; float64x2_t im; };
    complex128x2_t() : r{}, i{} {}
    complex128x2_t(float64x2_t re_, float64x2_t im_) : r(re_), i(im_) {}
};

/* ------------------------------------------------------------------ */
/* MF_IS_* macros injected by c_migrator for C scalar comparisons.    */
/* In C++ they just use the native == operator on MultiFloat.         */
/* ------------------------------------------------------------------ */

#define MF_IS_ZERO(x) ((x) == float64x2_t{0.0})
#define MF_IS_ONE(x)  ((x) == float64x2_t{1.0})

/* ------------------------------------------------------------------ */
/* Mixed-type operators.                                              */
/*                                                                    */
/* The upstream class defines only same-type operators plus a         */
/* non-explicit constructor from ``double``, so ``x op 1.0`` works    */
/* via implicit conversion. The reversed form ``1.0 op x`` needs      */
/* these free functions — PBLAS uses idioms like ``ONE + ssq*...``    */
/* with ONE=1.0 (double).                                             */
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
