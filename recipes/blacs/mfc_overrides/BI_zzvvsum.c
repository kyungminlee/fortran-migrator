#include "Bdef.h"
#include "multifloats_bridge.h"
void BI_zzvvsum(Int N, char *vec1, char *vec2)
{
   auto *v1 = (float64x2_t*)vec1;
   auto *v2 = (float64x2_t*)vec2;
   N *= 2;
   for (Int k = 0; k < N; k++) v1[k] = v1[k] + v2[k];
}
