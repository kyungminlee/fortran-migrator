/* BI_zzvvamn.c -- multifloats replacement for BI_zvvamn.
 *
 * Mirror of BI_zzvvamx; update when |v2| < |v1|.
 */
#include "Bdef.h"
#include "multifloats_c.h"

void BI_zzvvamn(Int N, char *vec1, char *vec2)
{
   complex128x2_t *v1 = (complex128x2_t *) vec1;
   complex128x2_t *v2 = (complex128x2_t *) vec2;
   BI_DistType *dist1, *dist2;
   Int i, k;

   k = N * sizeof(complex128x2_t);
   i = k % sizeof(BI_DistType);
   if (i) k += sizeof(BI_DistType) - i;
   dist1 = (BI_DistType *) &vec1[k];
   dist2 = (BI_DistType *) &vec2[k];

   for (k = 0; k < N; k++)
   {
      int c = mf_cmp(mf_cabs1(v1[k]), mf_cabs1(v2[k]));
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
