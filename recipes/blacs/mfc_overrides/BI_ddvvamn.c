#include "Bdef.h"
#include "multifloats_bridge.h"
void BI_ddvvamn(Int N, char *vec1, char *vec2)
{
   auto *v1 = (float64x2_t*)vec1;
   auto *v2 = (float64x2_t*)vec2;
   BI_DistType *dist1, *dist2;
   Int i, k;
   k = N * sizeof(float64x2_t);
   i = k % sizeof(BI_DistType);
   if (i) k += sizeof(BI_DistType) - i;
   dist1 = (BI_DistType *) &vec1[k];
   dist2 = (BI_DistType *) &vec2[k];
   for (k = 0; k < N; k++) {
      if (mf_abs(v1[k]) > mf_abs(v2[k])) {
         v1[k] = v2[k]; dist1[k] = dist2[k];
      } else if (mf_abs(v1[k]) == mf_abs(v2[k])) {
         if (dist1[k] > dist2[k]) { v1[k] = v2[k]; dist1[k] = dist2[k]; }
      }
   }
}
