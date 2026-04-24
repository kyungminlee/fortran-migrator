! ddlamch shim for the multifloats tests — see
! tests/lapack/target_kind16/lamch_shim.f90 for context.
!
! Unlike kind10/kind16 variants, the multifloats scalar type
! (float64x2) does not support Fortran's EPSILON/TINY/HUGE intrinsics
! directly, so the machine constants are built from the well-known
! double-double arithmetic identities:
!   eps_dd = 2^-104     ≈ 4.93e-32
!   tiny_dd ≈ tiny(1.0_dp) * 2^53   (denormals-free range)
!   huge_dd ≈ huge(1.0_dp)

function ddlamch(cmach) result(r)
    use multifloats, only: float64x2
    implicit none
    character(len=1), intent(in) :: cmach
    type(float64x2) :: r
    real(kind=8), parameter :: dzero = 0.0_8, done = 1.0_8
    real(kind=8) :: eps_dd_hi, sfmin_hi

    eps_dd_hi = 2.0_8**(-104)
    sfmin_hi  = tiny(dzero) * 2.0_8**53

    select case (cmach)
    case ('E', 'e'); r%limbs(1) = eps_dd_hi;         r%limbs(2) = dzero
    case ('S', 's'); r%limbs(1) = sfmin_hi;          r%limbs(2) = dzero
    case ('B', 'b'); r%limbs(1) = 2.0_8;             r%limbs(2) = dzero
    case ('P', 'p'); r%limbs(1) = eps_dd_hi * 2.0_8; r%limbs(2) = dzero
    case ('N', 'n'); r%limbs(1) = 106.0_8;           r%limbs(2) = dzero
    case ('R', 'r'); r%limbs(1) = done;              r%limbs(2) = dzero
    case ('M', 'm'); r%limbs(1) = real(minexponent(dzero) + 53, 8); r%limbs(2) = dzero
    case ('U', 'u'); r%limbs(1) = sfmin_hi;          r%limbs(2) = dzero
    case ('L', 'l'); r%limbs(1) = real(maxexponent(dzero),      8); r%limbs(2) = dzero
    case ('O', 'o'); r%limbs(1) = huge(dzero);       r%limbs(2) = dzero
    case default;    r%limbs(1) = dzero;             r%limbs(2) = dzero
    end select
end function ddlamch

function ddroundup_lwork(lwork) result(r)
    use multifloats, only: float64x2
    implicit none
    integer, intent(in) :: lwork
    type(float64x2) :: r
    r%limbs(1) = real(lwork, kind=8)
    r%limbs(2) = 0.0_8
    if (int(r%limbs(1)) < lwork) then
        r%limbs(1) = r%limbs(1) * (1.0_8 + 2.0_8**(-52))
    end if
end function ddroundup_lwork
