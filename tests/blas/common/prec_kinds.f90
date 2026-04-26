module prec_kinds
    implicit none
    private
    public :: ep, dp

    ! Universal reference precision for the differential precision
    ! tests — every comparison happens in REAL(KIND=ep) regardless of
    ! which target is being tested.
    integer, parameter :: ep = 16

    ! Double precision — used inside multifloats wrappers to bridge
    ! REAL(KIND=8) and TYPE(real64x2).
    integer, parameter :: dp = 8
end module prec_kinds
