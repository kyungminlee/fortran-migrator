/* pb_cttypeset.c -- multifloats replacement for PB_Cdtypeset.
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

PBTYP_T * PB_Cttypeset(void)
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

   TypeStruct.Cgesd2d   = Ctgesd2d;
   TypeStruct.Cgerv2d   = Ctgerv2d;
   TypeStruct.Cgebs2d   = Ctgebs2d;
   TypeStruct.Cgebr2d   = Ctgebr2d;
   TypeStruct.Cgsum2d   = Ctgsum2d;

   TypeStruct.Fmmadd    = tmmadd_;
   TypeStruct.Fmmcadd   = tmmcadd_;
   TypeStruct.Fmmtadd   = tmmtadd_;
   TypeStruct.Fmmtcadd  = tmmtcadd_;
   TypeStruct.Fmmdda    = tmmdda_;
   TypeStruct.Fmmddac   = tmmddac_;
   TypeStruct.Fmmddat   = tmmddat_;
   TypeStruct.Fmmddact  = tmmddact_;

   TypeStruct.Fcshft    = tcshft_;
   TypeStruct.Frshft    = trshft_;

   TypeStruct.Fvvdotu   = tvvdot_;
   TypeStruct.Fvvdotc   = tvvdot_;

   TypeStruct.Fset      = tset_;

   TypeStruct.Ftzpad    = ttzpad_;
   TypeStruct.Ftzpadcpy = ttzpadcpy_;
   TypeStruct.Ftzscal   = ttzscal_;
   TypeStruct.Fhescal   = ttzscal_;
   TypeStruct.Ftzcnjg   = ttzscal_;

   TypeStruct.Faxpy     = taxpy_;
   TypeStruct.Fcopy     = tcopy_;
   TypeStruct.Fswap     = tswap_;

   TypeStruct.Fgemv     = tgemv_;
   TypeStruct.Fsymv     = tsymv_;
   TypeStruct.Fhemv     = tsymv_;
   TypeStruct.Ftrmv     = ttrmv_;
   TypeStruct.Ftrsv     = ttrsv_;
   TypeStruct.Fagemv    = tagemv_;
   TypeStruct.Fasymv    = tasymv_;
   TypeStruct.Fahemv    = tasymv_;
   TypeStruct.Fatrmv    = tatrmv_;

   TypeStruct.Fgerc     = tger_;
   TypeStruct.Fgeru     = tger_;
   TypeStruct.Fsyr      = tsyr_;
   TypeStruct.Fher      = tsyr_;
   TypeStruct.Fsyr2     = tsyr2_;
   TypeStruct.Fher2     = tsyr2_;

   TypeStruct.Fgemm     = tgemm_;
   TypeStruct.Fsymm     = tsymm_;
   TypeStruct.Fhemm     = tsymm_;
   TypeStruct.Fsyrk     = tsyrk_;
   TypeStruct.Fherk     = tsyrk_;
   TypeStruct.Fsyr2k    = tsyr2k_;
   TypeStruct.Fher2k    = tsyr2k_;
   TypeStruct.Ftrmm     = ttrmm_;
   TypeStruct.Ftrsm     = ttrsm_;

   return( &TypeStruct );
}
