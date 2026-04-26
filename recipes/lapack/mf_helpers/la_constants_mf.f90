!> \brief \b LA_CONSTANTS_MF defines scaling constants for the multifloats
!> double-double (real64x2) precision, complementing LA_CONSTANTS.
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
   use multifloats, only: real64x2, cmplx64x2, &
                          DD_ZERO, DD_HALF, DD_ONE, DD_TWO, DD_EIGHT, &
                          DD_SAFMIN, DD_SAFMAX, &
                          DD_TSML, DD_TBIG, DD_SSML, DD_SBIG, &
                          DD_RTMIN, DD_RTMAX
   implicit none
   private
   public :: real64x2, cmplx64x2

! =====================================================================
!  Multifloats double-double (~106-bit) constants
! =====================================================================

!  Standard real constants
   type(real64x2), parameter, public :: ddzero  = DD_ZERO
   type(real64x2), parameter, public :: ddhalf  = DD_HALF
   type(real64x2), parameter, public :: ddone   = DD_ONE
   type(real64x2), parameter, public :: ddtwo   = DD_TWO
   ! Use named-component structure constructor (``limbs=...``) so that
   ! the compiler binds these initializers to the structure constructor
   ! of real64x2 rather than to the overloaded ``real64x2(...)``
   ! generic interface — the latter is a function call and is therefore
   ! illegal in a PARAMETER initializer.
   type(real64x2), parameter, public :: ddthree = real64x2(limbs=[3.0d0, 0.0d0])
   type(real64x2), parameter, public :: ddfour  = real64x2(limbs=[4.0d0, 0.0d0])
   type(real64x2), parameter, public :: ddeight = DD_EIGHT
   type(real64x2), parameter, public :: ddten   = real64x2(limbs=[10.0d0, 0.0d0])

!  Complex constants. Must use named-component structure constructor
!  syntax (``re=`` / ``im=``) so the compiler picks the structure
!  constructor and not the overloaded ``cmplx64x2`` interface
!  procedures (which are not allowed in PARAMETER initializers).
   type(cmplx64x2), parameter, public :: zzzero = &
      cmplx64x2(re=DD_ZERO, im=DD_ZERO)
   type(cmplx64x2), parameter, public :: zzhalf = &
      cmplx64x2(re=DD_HALF, im=DD_ZERO)
   type(cmplx64x2), parameter, public :: zzone  = &
      cmplx64x2(re=DD_ONE,  im=DD_ZERO)

   character*2, parameter, public :: ddprefix = 'DD'
   character*2, parameter, public :: zzprefix = 'ZZ'

!  Scaling constants (mirror la_constants.f90 names with dd prefix)
   type(real64x2), parameter, public :: ddsafmin = DD_SAFMIN
   type(real64x2), parameter, public :: ddsafmax = DD_SAFMAX
   type(real64x2), parameter, public :: ddrtmin  = DD_RTMIN
   type(real64x2), parameter, public :: ddrtmax  = DD_RTMAX

!  Blue's scaling constants
   type(real64x2), parameter, public :: ddtsml = DD_TSML
   type(real64x2), parameter, public :: ddtbig = DD_TBIG
   type(real64x2), parameter, public :: ddssml = DD_SSML
   type(real64x2), parameter, public :: ddsbig = DD_SBIG

end module LA_CONSTANTS_MF
