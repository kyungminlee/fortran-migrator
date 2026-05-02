/*
 * Target-aware real-arithmetic shim for the C-side MUMPS tests.
 *
 * The C tests (test_*_c_basic.c, test_*_c_sym.c) need a real type wide
 * enough to interoperate with the migrated bridge's typedef of
 * DMUMPS_REAL.  Each precision target ships its own
 * tests/mumps/target_<TARGET_NAME>/c/include/mumps_c_types.h that
 * defines DMUMPS_REAL appropriately:
 *
 *   kind16       — __float128
 *   kind10       — long double
 *   multifloats  — POD ``mumps_float64x2`` (no scalar arithmetic in C)
 *
 * This header picks the matching test-side type plus its literal suffix
 * macro, math helpers, epsilon, and snprintf wrapper.  Selected via
 * the ``TEST_TARGET_<NAME>`` define passed by the cmake glue.
 *
 * Multifloats is intentionally not supported here — its real type is a
 * POD struct with no operator overloading in plain C.  The cmake glue
 * skips C tests for that target until a portable harness lands.
 */

#ifndef TEST_REAL_COMPAT_H
#define TEST_REAL_COMPAT_H

#include <stdio.h>

#if defined(TEST_TARGET_KIND16)
#  include <quadmath.h>
   typedef __float128 test_real;
#  define TR_LIT(x)      x##q
#  define TR_FABS(x)     fabsq(x)
#  define TR_SQRT(x)     sqrtq(x)
#  define TR_EPS         FLT128_EPSILON
#  define TR_MIN         FLT128_MIN
   static inline int test_real_snprintf(char *buf, size_t n, test_real x) {
       return quadmath_snprintf(buf, n, "%.6Qe", x);
   }
#elif defined(TEST_TARGET_KIND10)
#  include <float.h>
#  include <math.h>
   typedef long double test_real;
#  define TR_LIT(x)      x##L
#  define TR_FABS(x)     fabsl(x)
#  define TR_SQRT(x)     sqrtl(x)
#  define TR_EPS         LDBL_EPSILON
#  define TR_MIN         LDBL_MIN
   static inline int test_real_snprintf(char *buf, size_t n, test_real x) {
       return snprintf(buf, n, "%.6Le", x);
   }
#else
#  error "test_real_compat.h: define TEST_TARGET_KIND16 or TEST_TARGET_KIND10"
#endif

#ifndef TEST_TARGET_NAME
#  error "test_real_compat.h: define TEST_TARGET_NAME (e.g. \"kind16\")"
#endif

#endif /* TEST_REAL_COMPAT_H */
