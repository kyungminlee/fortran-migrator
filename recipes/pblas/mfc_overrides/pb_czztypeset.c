/* pb_czztypeset.c -- multifloats replacement for PB_Cztypeset.
 *
 * Mirror of pb_cddtypeset.c for the complex double-double variant.
 * Initializes the static zero/one/negone constants in array form
 * (cmplxDD = float64x2_t[2]) and points the function pointers at
 * the migrated complex BLACS / BLAS / PBBLAS / PTZBLAS routines.
 */
#include "pblas.h"
#include "PBpblas.h"
#include "PBtools.h"
#include "PBblacs.h"
#include "PBblas.h"
#include "multifloats_c.h"

PBTYP_T * PB_Czztypeset(void)
{
   static Int     setup = 0;
   static PBTYP_T TypeStruct;
   static cmplxDD zero, one, negone;

   if( setup ) return( &TypeStruct );
   setup = 1;

   TypeStruct.type = DCPLX;
   TypeStruct.usiz = sizeof(float64x2_t);
   TypeStruct.size = sizeof(cmplxDD);

   zero[REAL_PART].limbs[0]   = 0.0; zero[REAL_PART].limbs[1]   = 0.0;
   zero[IMAG_PART].limbs[0]   = 0.0; zero[IMAG_PART].limbs[1]   = 0.0;
   one[REAL_PART].limbs[0]    = 1.0; one[REAL_PART].limbs[1]    = 0.0;
   one[IMAG_PART].limbs[0]    = 0.0; one[IMAG_PART].limbs[1]    = 0.0;
   negone[REAL_PART].limbs[0] = -1.0; negone[REAL_PART].limbs[1] = 0.0;
   negone[IMAG_PART].limbs[0] = 0.0; negone[IMAG_PART].limbs[1] = 0.0;

   TypeStruct.zero      = ((char *) zero);
   TypeStruct.one       = ((char *) one);
   TypeStruct.negone    = ((char *) negone);

   TypeStruct.Cgesd2d   = Czzgesd2d;
   TypeStruct.Cgerv2d   = Czzgerv2d;
   TypeStruct.Cgebs2d   = Czzgebs2d;
   TypeStruct.Cgebr2d   = Czzgebr2d;
   TypeStruct.Cgsum2d   = Czzgsum2d;

   TypeStruct.Fmmadd    = zzmmadd_;
   TypeStruct.Fmmcadd   = zzmmcadd_;
   TypeStruct.Fmmtadd   = zzmmtadd_;
   TypeStruct.Fmmtcadd  = zzmmtcadd_;
   TypeStruct.Fmmdda    = zzmmdda_;
   TypeStruct.Fmmddac   = zzmmddac_;
   TypeStruct.Fmmddat   = zzmmddat_;
   TypeStruct.Fmmddact  = zzmmddact_;

   TypeStruct.Fcshft    = zzcshft_;
   TypeStruct.Frshft    = zzrshft_;

   TypeStruct.Fvvdotu   = zzvvdotu_;
   TypeStruct.Fvvdotc   = zzvvdotc_;

   TypeStruct.Fset      = zzset_;

   TypeStruct.Ftzpad    = zztzpad_;
   TypeStruct.Ftzpadcpy = zztzpadcpy_;
   TypeStruct.Ftzscal   = zztzscal_;
   TypeStruct.Fhescal   = zzhescal_;
   TypeStruct.Ftzcnjg   = zztzcnjg_;

   TypeStruct.Faxpy     = zzaxpy_;
   TypeStruct.Fcopy     = zzcopy_;
   TypeStruct.Fswap     = zzswap_;

   TypeStruct.Fgemv     = zzgemv_;
   TypeStruct.Fsymv     = zzsymv_;
   TypeStruct.Fhemv     = zzhemv_;
   TypeStruct.Ftrmv     = zztrmv_;
   TypeStruct.Ftrsv     = zztrsv_;
   TypeStruct.Fagemv    = zzagemv_;
   TypeStruct.Fasymv    = zzasymv_;
   TypeStruct.Fahemv    = zzahemv_;
   TypeStruct.Fatrmv    = zzatrmv_;

   TypeStruct.Fgerc     = zzgerc_;
   TypeStruct.Fgeru     = zzgeru_;
   TypeStruct.Fsyr      = zzsyr_;
   TypeStruct.Fher      = zzher_;
   TypeStruct.Fsyr2     = zzsyr2_;
   TypeStruct.Fher2     = zzher2_;

   TypeStruct.Fgemm     = zzgemm_;
   TypeStruct.Fsymm     = zzsymm_;
   TypeStruct.Fhemm     = zzhemm_;
   TypeStruct.Fsyrk     = zzsyrk_;
   TypeStruct.Fherk     = zzherk_;
   TypeStruct.Fsyr2k    = zzsyr2k_;
   TypeStruct.Fher2k    = zzher2k_;
   TypeStruct.Ftrmm     = zztrmm_;
   TypeStruct.Ftrsm     = zztrsm_;

   return( &TypeStruct );
}
