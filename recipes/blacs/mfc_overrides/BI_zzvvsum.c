/* BI_zzvvsum.c -- multifloats replacement for BI_zvvsum.
 *
 * BI_zvvsum treats the complex vector as a flat array of 2N doubles.
 * We follow the same trick: treat it as 2N float64x2_t (real + imag
 * limbs stored consecutively, matching the complex128x2_t layout of
 * (re, im) where each is a float64x2_t).
 */
#include "Bdef.h"
#include "multifloats_c.h"

void BI_zzvvsum(Int N, char *vec1, char *vec2)
{
   float64x2_t *v1 = (float64x2_t *) vec1;
   float64x2_t *v2 = (float64x2_t *) vec2;
   Int k;
   N *= 2;
   for (k = 0; k < N; k++) v1[k] = mf_add(v1[k], v2[k]);
}
