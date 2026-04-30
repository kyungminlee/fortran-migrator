module prec_kinds
    implicit none
    private
    public :: ep, dp

    ! Universal reference precision for the MUMPS differential precision
    ! tests — every comparison happens in REAL(KIND=ep) regardless of
    ! which migrator target's qmumps archive is being exercised.
    integer, parameter :: ep = 16

    ! Double precision — used for legacy interfaces that still need it
    ! (e.g. MPI handles which are INTEGER == C int regardless of the
    ! Fortran working precision).
    integer, parameter :: dp = 8
end module prec_kinds
