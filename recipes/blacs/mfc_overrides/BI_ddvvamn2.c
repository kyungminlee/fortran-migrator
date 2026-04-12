#include "Bdef.h"
#include "multifloats_bridge.h"
void BI_ddvvamn2(Int N, char *vec1, char *vec2)
{
   auto *v1 = (float64x2_t*)vec1;
   auto *v2 = (float64x2_t*)vec2;
   for (Int k = 0; k != N; k++) {
      if (mf_abs(v1[k]) > mf_abs(v2[k])) v1[k] = v2[k];
      else if (mf_abs(v1[k]) == mf_abs(v2[k]))
         if (v1[k] < v2[k]) v1[k] = v2[k];
   }
}
