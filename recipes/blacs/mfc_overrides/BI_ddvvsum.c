/* BI_ddvvsum.c -- multifloats replacement for BI_dvvsum.
 *
 * The migrated-by-regex version would emit
 *     v1[k] += v2[k];
 * which is wrong for double-double because C struct assignment-add is
 * not defined on struct types. This hand-written replacement uses the
 * mf_add primitive from libmfc, which implements Knuth two-sum plus
 * renormalization for correct dd addition.
 */
#include "Bdef.h"
#include "multifloats_c.h"

void BI_ddvvsum(Int N, char *vec1, char *vec2)
{
   float64x2_t *v1 = (float64x2_t *) vec1;
   float64x2_t *v2 = (float64x2_t *) vec2;
   Int k;
   for (k = 0; k < N; k++) v1[k] = mf_add(v1[k], v2[k]);
}
