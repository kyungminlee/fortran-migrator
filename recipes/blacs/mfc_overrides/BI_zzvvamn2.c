/* BI_zzvvamn2.c -- multifloats replacement for BI_zvvamn2.
 *
 * Mirror of BI_zzvvamx2; pick the smaller |z|.
 */
#include "Bdef.h"
#include "multifloats_c.h"

void BI_zzvvamn2(Int N, char *vec1, char *vec2)
{
   complex128x2_t *v1 = (complex128x2_t *) vec1;
   complex128x2_t *v2 = (complex128x2_t *) vec2;
   Int k;

   for (k = 0; k < N; k++)
   {
      int c = mf_cmp(mf_cabs1(v1[k]), mf_cabs1(v2[k]));
      if (c > 0) v1[k] = v2[k];
      else if (c == 0)
      {
         int cr = mf_cmp(v1[k].re, v2[k].re);
         if (cr != 0)
         {
            if (cr < 0) v1[k] = v2[k];
         }
         else
         {
            if (mf_cmp(v1[k].im, v2[k].im) < 0) v1[k] = v2[k];
         }
      }
   }
}
