#include "Bdef.h"
#include "multifloats_bridge.h"
void BI_zzvvamx2(Int N, char *vec1, char *vec2)
{
   auto *v1 = (complex128x2_t*)vec1;
   auto *v2 = (complex128x2_t*)vec2;
   for (Int k = 0; k < N; k++) {
      if (mf_cabs1(v1[k]) < mf_cabs1(v2[k])) v1[k] = v2[k];
      else if (mf_cabs1(v1[k]) == mf_cabs1(v2[k])) {
         if (v1[k].re != v2[k].re) {
            if (v1[k].re < v2[k].re) v1[k] = v2[k];
         } else {
            if (v1[k].im < v2[k].im) v1[k] = v2[k];
         }
      }
   }
}
