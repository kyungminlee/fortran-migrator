/* pb_cddtypeset.c -- multifloats replacement for PB_Cdtypeset.
 *
 * The migrated-by-regex version cannot initialize the static
 * zero/one/negone constants with C operators on the float64x2_t
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
#include "multifloats_c.h"

PBTYP_T * PB_Cddtypeset(void)
{
   static Int     setup = 0;
   static PBTYP_T TypeStruct;
   static float64x2_t zero, one, negone;

   if( setup ) return( &TypeStruct );
   setup = 1;

   TypeStruct.type = DREAL;
   TypeStruct.usiz = sizeof(float64x2_t);
   TypeStruct.size = sizeof(float64x2_t);

   zero.limbs[0]   = 0.0;  zero.limbs[1]   = 0.0;
   one.limbs[0]    = 1.0;  one.limbs[1]    = 0.0;
   negone.limbs[0] = -1.0; negone.limbs[1] = 0.0;

   TypeStruct.zero      = (char *) (&zero);
   TypeStruct.one       = (char *) (&one);
   TypeStruct.negone    = (char *) (&negone);

   TypeStruct.Cgesd2d   = Cddgesd2d;
   TypeStruct.Cgerv2d   = Cddgerv2d;
   TypeStruct.Cgebs2d   = Cddgebs2d;
   TypeStruct.Cgebr2d   = Cddgebr2d;
   TypeStruct.Cgsum2d   = Cddgsum2d;

   TypeStruct.Fmmadd    = ddmmadd_;
   TypeStruct.Fmmcadd   = ddmmcadd_;
   TypeStruct.Fmmtadd   = ddmmtadd_;
   TypeStruct.Fmmtcadd  = ddmmtcadd_;
   TypeStruct.Fmmdda    = ddmmdda_;
   TypeStruct.Fmmddac   = ddmmddac_;
   TypeStruct.Fmmddat   = ddmmddat_;
   TypeStruct.Fmmddact  = ddmmddact_;

   TypeStruct.Fcshft    = ddcshft_;
   TypeStruct.Frshft    = ddrshft_;

   TypeStruct.Fvvdotu   = ddvvdot_;
   TypeStruct.Fvvdotc   = ddvvdot_;

   TypeStruct.Fset      = ddset_;

   TypeStruct.Ftzpad    = ddtzpad_;
   TypeStruct.Ftzpadcpy = ddtzpadcpy_;
   TypeStruct.Ftzscal   = ddtzscal_;
   TypeStruct.Fhescal   = ddtzscal_;
   TypeStruct.Ftzcnjg   = ddtzscal_;

   TypeStruct.Faxpy     = ddaxpy_;
   TypeStruct.Fcopy     = ddcopy_;
   TypeStruct.Fswap     = ddswap_;

   TypeStruct.Fgemv     = ddgemv_;
   TypeStruct.Fsymv     = ddsymv_;
   TypeStruct.Fhemv     = ddsymv_;
   TypeStruct.Ftrmv     = ddtrmv_;
   TypeStruct.Ftrsv     = ddtrsv_;
   TypeStruct.Fagemv    = ddagemv_;
   TypeStruct.Fasymv    = ddasymv_;
   TypeStruct.Fahemv    = ddasymv_;
   TypeStruct.Fatrmv    = ddatrmv_;

   TypeStruct.Fgerc     = ddger_;
   TypeStruct.Fgeru     = ddger_;
   TypeStruct.Fsyr      = ddsyr_;
   TypeStruct.Fher      = ddsyr_;
   TypeStruct.Fsyr2     = ddsyr2_;
   TypeStruct.Fher2     = ddsyr2_;

   TypeStruct.Fgemm     = ddgemm_;
   TypeStruct.Fsymm     = ddsymm_;
   TypeStruct.Fhemm     = ddsymm_;
   TypeStruct.Fsyrk     = ddsyrk_;
   TypeStruct.Fherk     = ddsyrk_;
   TypeStruct.Fsyr2k    = ddsyr2k_;
   TypeStruct.Fher2k    = ddsyr2k_;
   TypeStruct.Ftrmm     = ddtrmm_;
   TypeStruct.Ftrsm     = ddtrsm_;

   return( &TypeStruct );
}
