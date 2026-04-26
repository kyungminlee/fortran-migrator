/* pb_cvtypeset.c -- multifloats replacement for PB_Cztypeset.
 *
 * Mirror of pb_cttypeset.c for the complex double-double variant.
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

PBTYP_T * PB_Cvtypeset(void)
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

   TypeStruct.Cgesd2d   = Cvgesd2d;
   TypeStruct.Cgerv2d   = Cvgerv2d;
   TypeStruct.Cgebs2d   = Cvgebs2d;
   TypeStruct.Cgebr2d   = Cvgebr2d;
   TypeStruct.Cgsum2d   = Cvgsum2d;

   TypeStruct.Fmmadd    = vmmadd_;
   TypeStruct.Fmmcadd   = vmmcadd_;
   TypeStruct.Fmmtadd   = vmmtadd_;
   TypeStruct.Fmmtcadd  = vmmtcadd_;
   TypeStruct.Fmmdda    = vmmdda_;
   TypeStruct.Fmmddac   = vmmddac_;
   TypeStruct.Fmmddat   = vmmddat_;
   TypeStruct.Fmmddact  = vmmddact_;

   TypeStruct.Fcshft    = vcshft_;
   TypeStruct.Frshft    = vrshft_;

   TypeStruct.Fvvdotu   = vvvdotu_;
   TypeStruct.Fvvdotc   = vvvdotc_;

   TypeStruct.Fset      = vset_;

   TypeStruct.Ftzpad    = vtzpad_;
   TypeStruct.Ftzpadcpy = vtzpadcpy_;
   TypeStruct.Ftzscal   = vtzscal_;
   TypeStruct.Fhescal   = vhescal_;
   TypeStruct.Ftzcnjg   = vtzcnjg_;

   TypeStruct.Faxpy     = vaxpy_;
   TypeStruct.Fcopy     = vcopy_;
   TypeStruct.Fswap     = vswap_;

   TypeStruct.Fgemv     = vgemv_;
   TypeStruct.Fsymv     = vsymv_;
   TypeStruct.Fhemv     = vhemv_;
   TypeStruct.Ftrmv     = vtrmv_;
   TypeStruct.Ftrsv     = vtrsv_;
   TypeStruct.Fagemv    = vagemv_;
   TypeStruct.Fasymv    = vasymv_;
   TypeStruct.Fahemv    = vahemv_;
   TypeStruct.Fatrmv    = vatrmv_;

   TypeStruct.Fgerc     = vgerc_;
   TypeStruct.Fgeru     = vgeru_;
   TypeStruct.Fsyr      = vsyr_;
   TypeStruct.Fher      = vher_;
   TypeStruct.Fsyr2     = vsyr2_;
   TypeStruct.Fher2     = vher2_;

   TypeStruct.Fgemm     = vgemm_;
   TypeStruct.Fsymm     = vsymm_;
   TypeStruct.Fhemm     = vhemm_;
   TypeStruct.Fsyrk     = vsyrk_;
   TypeStruct.Fherk     = vherk_;
   TypeStruct.Fsyr2k    = vsyr2k_;
   TypeStruct.Fher2k    = vher2k_;
   TypeStruct.Ftrmm     = vtrmm_;
   TypeStruct.Ftrsm     = vtrsm_;

   return( &TypeStruct );
}
