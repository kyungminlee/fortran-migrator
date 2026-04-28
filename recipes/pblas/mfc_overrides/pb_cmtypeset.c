/* pb_cmtypeset.c -- multifloats replacement for PB_Cdtypeset.
 *
 * The migrated-by-regex version cannot initialize the static
 * zero/one/negone constants with C operators on the float64x2
 * struct (e.g. ``zero = ZERO`` where ZERO is the macro 0.0). This
 * hand-written version uses compound literals from libmfc to set
 * the constants and points the function pointers at the migrated
 * BLACS / BLAS / PBBLAS / PTZBLAS routines that the multifloats
 * BLACS migration produced.
 */
#include "pblas.h"
#include "PBpblas.h"
#include "PBtools.h"
#include "PBblacs.h"
#include "PBblas.h"
#include "multifloats_bridge.h"

PBTYP_T * PB_Cmtypeset(void)
{
   static Int     setup = 0;
   static PBTYP_T TypeStruct;
   static float64x2 zero, one, negone;

   if( setup ) return( &TypeStruct );
   setup = 1;

   TypeStruct.type = DREAL;
   TypeStruct.usiz = sizeof(float64x2);
   TypeStruct.size = sizeof(float64x2);

   zero   = 0.0;
   one    = 1.0;
   negone = -1.0;

   TypeStruct.zero      = (char *) (&zero);
   TypeStruct.one       = (char *) (&one);
   TypeStruct.negone    = (char *) (&negone);

   TypeStruct.Cgesd2d   = Cmgesd2d;
   TypeStruct.Cgerv2d   = Cmgerv2d;
   TypeStruct.Cgebs2d   = Cmgebs2d;
   TypeStruct.Cgebr2d   = Cmgebr2d;
   TypeStruct.Cgsum2d   = Cmgsum2d;

   TypeStruct.Fmmadd    = mmmadd_;
   TypeStruct.Fmmcadd   = mmmcadd_;
   TypeStruct.Fmmtadd   = mmmtadd_;
   TypeStruct.Fmmtcadd  = mmmtcadd_;
   TypeStruct.Fmmdda    = mmmdda_;
   TypeStruct.Fmmddac   = mmmddac_;
   TypeStruct.Fmmddat   = mmmddat_;
   TypeStruct.Fmmddact  = mmmddact_;

   TypeStruct.Fcshft    = mcshft_;
   TypeStruct.Frshft    = mrshft_;

   TypeStruct.Fvvdotu   = mvvdot_;
   TypeStruct.Fvvdotc   = mvvdot_;

   TypeStruct.Fset      = mset_;

   TypeStruct.Ftzpad    = mtzpad_;
   TypeStruct.Ftzpadcpy = mtzpadcpy_;
   TypeStruct.Ftzscal   = mtzscal_;
   TypeStruct.Fhescal   = mtzscal_;
   TypeStruct.Ftzcnjg   = mtzscal_;

   TypeStruct.Faxpy     = maxpy_;
   TypeStruct.Fcopy     = mcopy_;
   TypeStruct.Fswap     = mswap_;

   TypeStruct.Fgemv     = mgemv_;
   TypeStruct.Fsymv     = msymv_;
   TypeStruct.Fhemv     = msymv_;
   TypeStruct.Ftrmv     = mtrmv_;
   TypeStruct.Ftrsv     = mtrsv_;
   TypeStruct.Fagemv    = magemv_;
   TypeStruct.Fasymv    = masymv_;
   TypeStruct.Fahemv    = masymv_;
   TypeStruct.Fatrmv    = matrmv_;

   TypeStruct.Fgerc     = mger_;
   TypeStruct.Fgeru     = mger_;
   TypeStruct.Fsyr      = msyr_;
   TypeStruct.Fher      = msyr_;
   TypeStruct.Fsyr2     = msyr2_;
   TypeStruct.Fher2     = msyr2_;

   TypeStruct.Fgemm     = mgemm_;
   TypeStruct.Fsymm     = msymm_;
   TypeStruct.Fhemm     = msymm_;
   TypeStruct.Fsyrk     = msyrk_;
   TypeStruct.Fherk     = msyrk_;
   TypeStruct.Fsyr2k    = msyr2k_;
   TypeStruct.Fher2k    = msyr2k_;
   TypeStruct.Ftrmm     = mtrmm_;
   TypeStruct.Ftrsm     = mtrsm_;

   return( &TypeStruct );
}
