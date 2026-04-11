/* BI_ddvvamx.c -- multifloats replacement for BI_dvvamx.
 *
 * The scalar ``diff = Rabs(v1) - Rabs(v2)'' sign test cannot be done
 * with C operators on a struct type, so we use mf_cmp(mf_abs(...)) to
 * compare the two double-double absolute values lexicographically.
 * Tie-breaking via the BI_DistType tail pointer is preserved verbatim.
 */
#include "Bdef.h"
#include "multifloats_c.h"

void BI_ddvvamx(Int N, char *vec1, char *vec2)
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
      if (c < 0)
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
