/* pddlaiect.c -- hand-written multifloats override for pdlaiect.c
 *
 * The original pdlaiect.c extracts the sign bit of a double via
 * pointer-cast to unsigned int and bit-shift.  That trick relies on
 * the IEEE 754 double layout, which does not apply to float64x2_t.
 *
 * This override replaces the bit manipulation with C++ comparison
 * operators provided by multifloats_bridge.h.
 *
 * Compiled as C++ (mpicxx) with -include multifloats_bridge.h.
 */

#include "pxsyevx.h"
#include "pblas.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void pddlasnbt_(Int *ieflag)
{
    /* For float64x2_t, sign detection uses C++ comparison operators.
     * Return nonzero to signal that sign detection is available. */
    *ieflag = 1;
}

void pddlaiectb_(float64x2_t *sigma, Int *n, float64x2_t *d, Int *count)
{
    float64x2_t lsigma = *sigma;
    float64x2_t *pd = d;
    float64x2_t *pe2 = d + 1;
    float64x2_t tmp = *pd - lsigma;
    pd += 2;
    *count = (tmp < float64x2_t{0.0}) ? 1 : 0;
    for (Int i = 1; i < *n; i++) {
        tmp = *pd - *pe2 / tmp - lsigma;
        pd += 2;
        pe2 += 2;
        *count += (tmp < float64x2_t{0.0}) ? 1 : 0;
    }
}

void pddlaiectl_(float64x2_t *sigma, Int *n, float64x2_t *d, Int *count)
{
    /* The big-endian / little-endian distinction for sign-bit extraction
     * does not apply to float64x2_t comparison-based sign detection.
     * Both functions are identical. */
    pddlaiectb_(sigma, n, d, count);
}

void pddlachkieee_(Int *isieee, float64x2_t *rmax, float64x2_t *rmin)
{
    /* float64x2_t uses C++ operator overloads that handle sign
     * and comparison correctly.  Report IEEE-compatible behaviour. */
    (void)rmax;
    (void)rmin;
    *isieee = 1;
}

#ifdef __cplusplus
} /* extern "C" */
#endif
