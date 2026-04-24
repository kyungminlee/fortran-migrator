! qlamch shim for the kind16 tests.
!
! The migrator (recipe/lapack.yaml) declares
! external/lapack-3.12.1/INSTALL as an `extra_symbol_dirs` entry, but
! `cmd_stage` does not copy and migrate those files — so DLAMCH (and
! its prefixed variants QLAMCH / ELAMCH / DDLAMCH) never appears in
! the migrated library. Without it, qlapack has unresolved references
! (from qlascl, qsyev, qgesvd, qgetrf2, …) and no LAPACK test links.
!
! This shim locally reconstructs QLAMCH using Fortran intrinsics at
! KIND=16 — the same structure as INSTALL/dlamch.f, just directly at
! the target precision. Remove this file once the migrator stages
! INSTALL/ properly.

real(kind=16) function qlamch(cmach) result(r)
    implicit none
    character(len=1), intent(in) :: cmach
    real(kind=16), parameter :: one = 1.0_16, zero = 0.0_16
    real(kind=16) :: eps, sfmin, small

    ! Assume rounding, not chopping. Always.
    eps = epsilon(zero) * 0.5_16

    select case (cmach)
    case ('E', 'e'); r = eps
    case ('S', 's')
        sfmin = tiny(zero)
        small = one / huge(zero)
        if (small >= sfmin) sfmin = small * (one + eps)
        r = sfmin
    case ('B', 'b'); r = real(radix(zero),        16)
    case ('P', 'p'); r = eps * real(radix(zero),  16)
    case ('N', 'n'); r = real(digits(zero),       16)
    case ('R', 'r'); r = one
    case ('M', 'm'); r = real(minexponent(zero),  16)
    case ('U', 'u'); r = tiny(zero)
    case ('L', 'l'); r = real(maxexponent(zero),  16)
    case ('O', 'o'); r = huge(zero)
    case default;    r = zero
    end select
end function qlamch

real(kind=16) function qroundup_lwork(lwork) result(r)
    implicit none
    integer, intent(in) :: lwork
    r = real(lwork, kind=16)
    if (int(r) < lwork) then
        r = r * (1.0_16 + epsilon(0.0_16))
    end if
end function qroundup_lwork
