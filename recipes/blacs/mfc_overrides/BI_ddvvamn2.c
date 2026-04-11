/* BI_ddvvamn2.c -- multifloats replacement for BI_dvvamn2.
 *
 * Tie-break by direct value compare, mirror of BI_ddvvamx2.
 */
#include "Bdef.h"
#include "multifloats_c.h"

void BI_ddvvamn2(Int N, char *vec1, char *vec2)
{
   float64x2_t *v1 = (float64x2_t *) vec1;
   float64x2_t *v2 = (float64x2_t *) vec2;
   Int k;

   for (k = 0; k != N; k++)
   {
      int c = mf_cmp(mf_abs(v1[k]), mf_abs(v2[k]));
      if (c > 0) v1[k] = v2[k];
      else if (c == 0)
      {
         if (mf_cmp(v1[k], v2[k]) < 0) v1[k] = v2[k];
      }
   }
}
