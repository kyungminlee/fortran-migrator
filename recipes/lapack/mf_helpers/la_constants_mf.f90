!> \brief \b LA_CONSTANTS_MF defines scaling constants for the multifloats
!> double-double (real64x2) precision, complementing LA_CONSTANTS.
!
!  =========== DOCUMENTATION ===========
!
!  Multifloats double-double (~106-bit mantissa) constants. Re-exports
!  named constants from the multifloats module under the migrator's
!  multifloats prefix (T for real, V for complex). The migrator
!  rewrites the USE-clause aliases in LAPACK source from
!  ``zero=>dzero`` (etc.) to ``zero=>tzero`` (etc.); the local LHS
!  alias is unchanged so the body of the routine compiles without
!  further substitution.
!
!  Naming conventions per migrator target:
!    KIND=8  (double):  D prefix (real), Z prefix (complex)
!    KIND=10 (extended): E prefix (real), Y prefix (complex)
!    KIND=16 (quad):    Q prefix (real), X prefix (complex)
!    multifloats:       T prefix (real), V prefix (complex)
!  Single-letter T/V are unused as routine prefixes anywhere in
!  upstream BLAS / LAPACK / ScaLAPACK, so the renamed namespace is
!  collision-free. (Earlier two-letter DD/ZZ prefixes collided with
!  ScaLAPACK's orphaned DDDOT wrapper, corrupting pddpotf2.)
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

!  Standard real constants (T-prefixed to match the migrator rename of
!  the corresponding D-prefixed source symbols).
   type(real64x2), parameter, public :: tzero  = DD_ZERO
   type(real64x2), parameter, public :: thalf  = DD_HALF
   type(real64x2), parameter, public :: tone   = DD_ONE
   type(real64x2), parameter, public :: ttwo   = DD_TWO
   ! Use named-component structure constructor (``limbs=...``) so that
   ! the compiler binds these initializers to the structure constructor
   ! of real64x2 rather than to the overloaded ``real64x2(...)``
   ! generic interface — the latter is a function call and is therefore
   ! illegal in a PARAMETER initializer.
   type(real64x2), parameter, public :: tthree = real64x2(limbs=[3.0d0, 0.0d0])
   type(real64x2), parameter, public :: tfour  = real64x2(limbs=[4.0d0, 0.0d0])
   type(real64x2), parameter, public :: teight = DD_EIGHT
   type(real64x2), parameter, public :: tten   = real64x2(limbs=[10.0d0, 0.0d0])

!  Complex constants (V-prefixed). Must use named-component structure
!  constructor syntax (``re=`` / ``im=``) so the compiler picks the
!  structure constructor and not the overloaded ``cmplx64x2`` interface
!  procedures (which are not allowed in PARAMETER initializers).
   type(cmplx64x2), parameter, public :: vzero = &
      cmplx64x2(re=DD_ZERO, im=DD_ZERO)
   type(cmplx64x2), parameter, public :: vhalf = &
      cmplx64x2(re=DD_HALF, im=DD_ZERO)
   type(cmplx64x2), parameter, public :: vone  = &
      cmplx64x2(re=DD_ONE,  im=DD_ZERO)

   character*1, parameter, public :: tprefix = 'T'
   character*1, parameter, public :: vprefix = 'V'

!  Scaling constants (mirror la_constants.f90 names with t prefix)
   type(real64x2), parameter, public :: tsafmin = DD_SAFMIN
   type(real64x2), parameter, public :: tsafmax = DD_SAFMAX
   type(real64x2), parameter, public :: trtmin  = DD_RTMIN
   type(real64x2), parameter, public :: trtmax  = DD_RTMAX

!  Blue's scaling constants
   type(real64x2), parameter, public :: ttsml = DD_TSML
   type(real64x2), parameter, public :: ttbig = DD_TBIG
   type(real64x2), parameter, public :: tssml = DD_SSML
   type(real64x2), parameter, public :: tsbig = DD_SBIG

end module LA_CONSTANTS_MF
