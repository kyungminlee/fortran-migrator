! elamch shim for the kind10 tests — see tests/lapack/target_kind16/lamch_shim.f90
! for context. Same structure as INSTALL/dlamch.f at KIND=10.

real(kind=10) function elamch(cmach) result(r)
    implicit none
    character(len=1), intent(in) :: cmach
    real(kind=10), parameter :: one = 1.0_10, zero = 0.0_10
    real(kind=10) :: eps, sfmin, small

    eps = epsilon(zero) * 0.5_10

    select case (cmach)
    case ('E', 'e'); r = eps
    case ('S', 's')
        sfmin = tiny(zero)
        small = one / huge(zero)
        if (small >= sfmin) sfmin = small * (one + eps)
        r = sfmin
    case ('B', 'b'); r = real(radix(zero),        10)
    case ('P', 'p'); r = eps * real(radix(zero),  10)
    case ('N', 'n'); r = real(digits(zero),       10)
    case ('R', 'r'); r = one
    case ('M', 'm'); r = real(minexponent(zero),  10)
    case ('U', 'u'); r = tiny(zero)
    case ('L', 'l'); r = real(maxexponent(zero),  10)
    case ('O', 'o'); r = huge(zero)
    case default;    r = zero
    end select
end function elamch

real(kind=10) function eroundup_lwork(lwork) result(r)
    implicit none
    integer, intent(in) :: lwork
    r = real(lwork, kind=10)
    if (int(r) < lwork) then
        r = r * (1.0_10 + epsilon(0.0_10))
    end if
end function eroundup_lwork
