/* pb_cwtypeset.c -- multifloats replacement for PB_Cztypeset.
 *
 * Mirror of pb_cmtypeset.c for the complex double-double variant.
 * Initializes the static zero/one/negone constants in array form
 * (cmplxDD = float64x2[2]) and points the function pointers at
 * the migrated complex BLACS / BLAS / PBBLAS / PTZBLAS routines.
 */
#include "pblas.h"
#include "PBpblas.h"
#include "PBtools.h"
#include "PBblacs.h"
#include "PBblas.h"
#include "multifloats_bridge.h"

PBTYP_T * PB_Cwtypeset(void)
{
   static Int     setup = 0;
   static PBTYP_T TypeStruct;
   static cmplxDD zero, one, negone;

   if( setup ) return( &TypeStruct );
   setup = 1;

   TypeStruct.type = DCPLX;
   TypeStruct.usiz = sizeof(float64x2);
   TypeStruct.size = sizeof(cmplxDD);

   zero[REAL_PART]   = 0.0;
   zero[IMAG_PART]   = 0.0;
   one[REAL_PART]    = 1.0;
   one[IMAG_PART]    = 0.0;
   negone[REAL_PART] = -1.0;
   negone[IMAG_PART] = 0.0;

   TypeStruct.zero      = ((char *) zero);
   TypeStruct.one       = ((char *) one);
   TypeStruct.negone    = ((char *) negone);

   TypeStruct.Cgesd2d   = Cwgesd2d;
   TypeStruct.Cgerv2d   = Cwgerv2d;
   TypeStruct.Cgebs2d   = Cwgebs2d;
   TypeStruct.Cgebr2d   = Cwgebr2d;
   TypeStruct.Cgsum2d   = Cwgsum2d;

   TypeStruct.Fmmadd    = wmmadd_;
   TypeStruct.Fmmcadd   = wmmcadd_;
   TypeStruct.Fmmtadd   = wmmtadd_;
   TypeStruct.Fmmtcadd  = wmmtcadd_;
   TypeStruct.Fmmdda    = wmmdda_;
   TypeStruct.Fmmddac   = wmmddac_;
   TypeStruct.Fmmddat   = wmmddat_;
   TypeStruct.Fmmddact  = wmmddact_;

   TypeStruct.Fcshft    = wcshft_;
   TypeStruct.Frshft    = wrshft_;

   TypeStruct.Fvvdotu   = wvvdotu_;
   TypeStruct.Fvvdotc   = wvvdotc_;

   TypeStruct.Fset      = wset_;

   TypeStruct.Ftzpad    = wtzpad_;
   TypeStruct.Ftzpadcpy = wtzpadcpy_;
   TypeStruct.Ftzscal   = wtzscal_;
   TypeStruct.Fhescal   = whescal_;
   TypeStruct.Ftzcnjg   = wtzcnjg_;

   TypeStruct.Faxpy     = waxpy_;
   TypeStruct.Fcopy     = wcopy_;
   TypeStruct.Fswap     = wswap_;

   TypeStruct.Fgemv     = wgemv_;
   TypeStruct.Fsymv     = wsymv_;
   TypeStruct.Fhemv     = whemv_;
   TypeStruct.Ftrmv     = wtrmv_;
   TypeStruct.Ftrsv     = wtrsv_;
   TypeStruct.Fagemv    = wagemv_;
   TypeStruct.Fasymv    = wasymv_;
   TypeStruct.Fahemv    = wahemv_;
   TypeStruct.Fatrmv    = watrmv_;

   TypeStruct.Fgerc     = wgerc_;
   TypeStruct.Fgeru     = wgeru_;
   TypeStruct.Fsyr      = wsyr_;
   TypeStruct.Fher      = wher_;
   TypeStruct.Fsyr2     = wsyr2_;
   TypeStruct.Fher2     = wher2_;

   TypeStruct.Fgemm     = wgemm_;
   TypeStruct.Fsymm     = wsymm_;
   TypeStruct.Fhemm     = whemm_;
   TypeStruct.Fsyrk     = wsyrk_;
   TypeStruct.Fherk     = wherk_;
   TypeStruct.Fsyr2k    = wsyr2k_;
   TypeStruct.Fher2k    = wher2k_;
   TypeStruct.Ftrmm     = wtrmm_;
   TypeStruct.Ftrsm     = wtrsm_;

   return( &TypeStruct );
}
