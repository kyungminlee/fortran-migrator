!> \brief \b LA_CONSTANTS_MF defines scaling constants for the multifloats
!> double-double (float64x2) precision, complementing LA_CONSTANTS.
!
!  =========== DOCUMENTATION ===========
!
!  Multifloats double-double (~106-bit mantissa) constants. Re-exports
!  named constants from the multifloats module under the w-prefixed
!  (real) and u-prefixed (complex) names. The migrator rewrites the
!  USE-clause aliases in LAPACK source from ``zero=>dzero`` (etc.) to
!  ``zero=>wzero`` (etc.); the local LHS alias is unchanged so the
!  body of the routine compiles without further substitution.
!
!  Naming conventions for multifloats are distinct from existing
!  precisions:
!    KIND=8  (double):  D prefix (real), Z prefix (complex)
!    KIND=10 (extended): E prefix (real), Y prefix (complex)
!    KIND=16 (quad):    Q prefix (real), X prefix (complex)
!    multifloats:       W prefix (real), U prefix (complex)
!  W stands for "wide" — wider precision than double via a different
!  mechanism than hardware quad.
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
   type(float64x2), parameter, public :: wzero  = MF_ZERO
   type(float64x2), parameter, public :: whalf  = MF_HALF
   type(float64x2), parameter, public :: wone   = MF_ONE
   type(float64x2), parameter, public :: wtwo   = MF_TWO
   ! Use named-component structure constructor (``limbs=...``) so that
   ! the compiler binds these initializers to the structure constructor
   ! of float64x2 rather than to the overloaded ``float64x2(...)``
   ! generic interface — the latter is a function call and is therefore
   ! illegal in a PARAMETER initializer.
   type(float64x2), parameter, public :: wthree = float64x2(limbs=[3.0d0, 0.0d0])
   type(float64x2), parameter, public :: wfour  = float64x2(limbs=[4.0d0, 0.0d0])
   type(float64x2), parameter, public :: weight = MF_EIGHT
   type(float64x2), parameter, public :: wten   = float64x2(limbs=[10.0d0, 0.0d0])

!  Complex constants. Must use named-component structure constructor
!  syntax (``re=`` / ``im=``) so the compiler picks the structure
!  constructor and not the overloaded ``complex128x2`` interface
!  procedures (which are not allowed in PARAMETER initializers).
   type(complex128x2), parameter, public :: uzero = &
      complex128x2(re=MF_ZERO, im=MF_ZERO)
   type(complex128x2), parameter, public :: uhalf = &
      complex128x2(re=MF_HALF, im=MF_ZERO)
   type(complex128x2), parameter, public :: uone  = &
      complex128x2(re=MF_ONE,  im=MF_ZERO)

   character*1, parameter, public :: wprefix = 'W'
   character*1, parameter, public :: uprefix = 'U'

!  Scaling constants (mirror la_constants.f90 names with w prefix)
   type(float64x2), parameter, public :: wsafmin = MF_SAFMIN
   type(float64x2), parameter, public :: wsafmax = MF_SAFMAX
   type(float64x2), parameter, public :: wrtmin  = MF_RTMIN
   type(float64x2), parameter, public :: wrtmax  = MF_RTMAX

!  Blue's scaling constants
   type(float64x2), parameter, public :: wtsml = MF_TSML
   type(float64x2), parameter, public :: wtbig = MF_TBIG
   type(float64x2), parameter, public :: wssml = MF_SSML
   type(float64x2), parameter, public :: wsbig = MF_SBIG

end module LA_CONSTANTS_MF
