/* BI_ddvvamn.c -- multifloats replacement for BI_dvvamn.
 *
 * Mirror of BI_ddvvamx with the sense of the abs compare flipped:
 * update when |v2| < |v1|. Tie-breaking via BI_DistType tail preserved.
 */
#include "Bdef.h"
#include "multifloats_c.h"

void BI_ddvvamn(Int N, char *vec1, char *vec2)
{
   float64x2_t *v1 = (float64x2_t *) vec1;
   float64x2_t *v2 = (float64x2_t *) vec2;
   BI_DistType *dist1, *dist2;
   Int i, k;

   k = N * sizeof(float64x2_t);
   i = k % sizeof(BI_DistType);
   if (i) k += sizeof(BI_DistType) - i;
   dist1 = (BI_DistType *) &vec1[k];
   dist2 = (BI_DistType *) &vec2[k];

   for (k = 0; k < N; k++)
   {
      int c = mf_cmp(mf_abs(v1[k]), mf_abs(v2[k]));
      if (c > 0)
      {
         v1[k] = v2[k];
         dist1[k] = dist2[k];
      }
      else if (c == 0)
      {
         if (dist1[k] > dist2[k])
         {
            v1[k] = v2[k];
            dist1[k] = dist2[k];
         }
      }
   }
}
