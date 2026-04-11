!> \brief \b LA_CONSTANTS_MF defines scaling constants for the multifloats
!> double-double (float64x2) precision, complementing LA_CONSTANTS.
!
!  =========== DOCUMENTATION ===========
!
!  Multifloats double-double (~106-bit mantissa) constants. Re-exports
!  named constants from the multifloats module under the dd-prefixed
!  (real) and zz-prefixed (complex) names. The migrator rewrites the
!  USE-clause aliases in LAPACK source from ``zero=>dzero`` (etc.) to
!  ``zero=>ddzero`` (etc.); the local LHS alias is unchanged so the
!  body of the routine compiles without further substitution.
!
!  Naming conventions for multifloats are distinct from existing
!  precisions:
!    KIND=8  (double):  D prefix (real), Z prefix (complex)
!    KIND=10 (extended): E prefix (real), Y prefix (complex)
!    KIND=16 (quad):    Q prefix (real), X prefix (complex)
!    multifloats:       DD prefix (real), ZZ prefix (complex)
!  Two-letter prefixes avoid the single-letter collisions that bit
!  the earlier ``W``/``U`` convention — LAPACK reserves ``W`` as the
!  workspace-size idiom (e.g. WLALSD in DGELSD), so renaming DLALSD
!  to WLALSD created a conflict with the existing local integer.
!
module LA_CONSTANTS_MF
   use multifloats, only: float64x2, complex128x2, &
                          MF_ZERO, MF_HALF, MF_ONE, MF_TWO, MF_EIGHT, &
                          MF_SAFMIN, MF_SAFMAX, &
                          MF_TSML, MF_TBIG, MF_SSML, MF_SBIG, &
                          MF_RTMIN, MF_RTMAX
   implicit none
   private
   public :: float64x2, complex128x2

! =====================================================================
!  Multifloats double-double (~106-bit) constants
! =====================================================================

!  Standard real constants
   type(float64x2), parameter, public :: ddzero  = MF_ZERO
   type(float64x2), parameter, public :: ddhalf  = MF_HALF
   type(float64x2), parameter, public :: ddone   = MF_ONE
   type(float64x2), parameter, public :: ddtwo   = MF_TWO
   ! Use named-component structure constructor (``limbs=...``) so that
   ! the compiler binds these initializers to the structure constructor
   ! of float64x2 rather than to the overloaded ``float64x2(...)``
   ! generic interface — the latter is a function call and is therefore
   ! illegal in a PARAMETER initializer.
   type(float64x2), parameter, public :: ddthree = float64x2(limbs=[3.0d0, 0.0d0])
   type(float64x2), parameter, public :: ddfour  = float64x2(limbs=[4.0d0, 0.0d0])
   type(float64x2), parameter, public :: ddeight = MF_EIGHT
   type(float64x2), parameter, public :: ddten   = float64x2(limbs=[10.0d0, 0.0d0])

!  Complex constants. Must use named-component structure constructor
!  syntax (``re=`` / ``im=``) so the compiler picks the structure
!  constructor and not the overloaded ``complex128x2`` interface
!  procedures (which are not allowed in PARAMETER initializers).
   type(complex128x2), parameter, public :: zzzero = &
      complex128x2(re=MF_ZERO, im=MF_ZERO)
   type(complex128x2), parameter, public :: zzhalf = &
      complex128x2(re=MF_HALF, im=MF_ZERO)
   type(complex128x2), parameter, public :: zzone  = &
      complex128x2(re=MF_ONE,  im=MF_ZERO)

   character*2, parameter, public :: ddprefix = 'DD'
   character*2, parameter, public :: zzprefix = 'ZZ'

!  Scaling constants (mirror la_constants.f90 names with dd prefix)
   type(float64x2), parameter, public :: ddsafmin = MF_SAFMIN
   type(float64x2), parameter, public :: ddsafmax = MF_SAFMAX
   type(float64x2), parameter, public :: ddrtmin  = MF_RTMIN
   type(float64x2), parameter, public :: ddrtmax  = MF_RTMAX

!  Blue's scaling constants
   type(float64x2), parameter, public :: ddtsml = MF_TSML
   type(float64x2), parameter, public :: ddtbig = MF_TBIG
   type(float64x2), parameter, public :: ddssml = MF_SSML
   type(float64x2), parameter, public :: ddsbig = MF_SBIG

end module LA_CONSTANTS_MF
