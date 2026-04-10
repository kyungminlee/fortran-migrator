






module multifloats
    implicit none
    private

    integer, parameter :: qp = 16  ! Quad precision
    integer, parameter :: dp = 8   ! Double precision
    integer, parameter :: sp = 4   ! Single precision

    type, public :: float64x2
        real(dp) :: limbs(2)
    end type float64x2

    type, public :: complex128x2
        type(float64x2) :: re, im
    end type complex128x2

    ! ================================================================
    ! Constructors
    ! ================================================================

    interface float64x2
        module procedure mf_from_dp
        module procedure mf_from_sp
        module procedure mf_from_int
        module procedure mf_from_char
        module procedure mf_from_cx
        module procedure mf_from_cdp
        module procedure mf_from_csp
    end interface float64x2

    interface complex128x2
        module procedure cx_from_mf_2
        module procedure cx_from_dp_2
        module procedure cx_from_mf_1
        module procedure cx_from_dp_1
    end interface complex128x2

    ! ================================================================
    ! Arithmetic operators (+, -, *, /)
    ! ================================================================

    interface operator(+)
        module procedure mf_add_mf
        module procedure mf_add_dp
        module procedure mf_add_sp
        module procedure mf_add_int
        module procedure mf_add_cx
        module procedure mf_add_cdp
        module procedure mf_add_csp
        module procedure dp_add_mf
        module procedure dp_add_cx
        module procedure sp_add_mf
        module procedure sp_add_cx
        module procedure int_add_mf
        module procedure int_add_cx
        module procedure cx_add_mf
        module procedure cx_add_dp
        module procedure cx_add_sp
        module procedure cx_add_int
        module procedure cx_add_cx
        module procedure cx_add_cdp
        module procedure cx_add_csp
        module procedure cdp_add_mf
        module procedure cdp_add_cx
        module procedure csp_add_mf
        module procedure csp_add_cx
    end interface
    public :: operator(+)

    interface operator(-)
        module procedure mf_sub_mf
        module procedure mf_sub_dp
        module procedure mf_sub_sp
        module procedure mf_sub_int
        module procedure mf_sub_cx
        module procedure mf_sub_cdp
        module procedure mf_sub_csp
        module procedure dp_sub_mf
        module procedure dp_sub_cx
        module procedure sp_sub_mf
        module procedure sp_sub_cx
        module procedure int_sub_mf
        module procedure int_sub_cx
        module procedure cx_sub_mf
        module procedure cx_sub_dp
        module procedure cx_sub_sp
        module procedure cx_sub_int
        module procedure cx_sub_cx
        module procedure cx_sub_cdp
        module procedure cx_sub_csp
        module procedure cdp_sub_mf
        module procedure cdp_sub_cx
        module procedure csp_sub_mf
        module procedure csp_sub_cx
        module procedure mf_neg
        module procedure cx_neg
    end interface
    public :: operator(-)

    interface operator(*)
        module procedure mf_mul_mf
        module procedure mf_mul_dp
        module procedure mf_mul_sp
        module procedure mf_mul_int
        module procedure mf_mul_cx
        module procedure mf_mul_cdp
        module procedure mf_mul_csp
        module procedure dp_mul_mf
        module procedure dp_mul_cx
        module procedure sp_mul_mf
        module procedure sp_mul_cx
        module procedure int_mul_mf
        module procedure int_mul_cx
        module procedure cx_mul_mf
        module procedure cx_mul_dp
        module procedure cx_mul_sp
        module procedure cx_mul_int
        module procedure cx_mul_cx
        module procedure cx_mul_cdp
        module procedure cx_mul_csp
        module procedure cdp_mul_mf
        module procedure cdp_mul_cx
        module procedure csp_mul_mf
        module procedure csp_mul_cx
    end interface
    public :: operator(*)

    interface operator(/)
        module procedure mf_div_mf
        module procedure mf_div_dp
        module procedure mf_div_sp
        module procedure mf_div_int
        module procedure mf_div_cx
        module procedure mf_div_cdp
        module procedure mf_div_csp
        module procedure dp_div_mf
        module procedure dp_div_cx
        module procedure sp_div_mf
        module procedure sp_div_cx
        module procedure int_div_mf
        module procedure int_div_cx
        module procedure cx_div_mf
        module procedure cx_div_dp
        module procedure cx_div_sp
        module procedure cx_div_int
        module procedure cx_div_cx
        module procedure cx_div_cdp
        module procedure cx_div_csp
        module procedure cdp_div_mf
        module procedure cdp_div_cx
        module procedure csp_div_mf
        module procedure csp_div_cx
    end interface
    public :: operator(/)

    ! ================================================================
    ! Power operator (**)
    ! ================================================================

    interface operator(**)
        module procedure mf_pow_mf
        module procedure mf_pow_dp
        module procedure mf_pow_sp
        module procedure mf_pow_int
        module procedure mf_pow_cx
        module procedure mf_pow_cdp
        module procedure mf_pow_csp
        module procedure dp_pow_mf
        module procedure dp_pow_cx
        module procedure sp_pow_mf
        module procedure sp_pow_cx
        module procedure int_pow_mf
        module procedure int_pow_cx
        module procedure cx_pow_mf
        module procedure cx_pow_dp
        module procedure cx_pow_sp
        module procedure cx_pow_int
        module procedure cx_pow_cx
        module procedure cx_pow_cdp
        module procedure cx_pow_csp
        module procedure cdp_pow_mf
        module procedure cdp_pow_cx
        module procedure csp_pow_mf
        module procedure csp_pow_cx
    end interface
    public :: operator(**)

    ! ================================================================
    ! Equality comparison (==, /=) — all type pairs
    ! ================================================================

    interface operator(==)
        module procedure mf_eq_mf
        module procedure mf_eq_dp
        module procedure mf_eq_sp
        module procedure mf_eq_int
        module procedure mf_eq_cx
        module procedure mf_eq_cdp
        module procedure mf_eq_csp
        module procedure dp_eq_mf
        module procedure dp_eq_cx
        module procedure sp_eq_mf
        module procedure sp_eq_cx
        module procedure int_eq_mf
        module procedure int_eq_cx
        module procedure cx_eq_mf
        module procedure cx_eq_dp
        module procedure cx_eq_sp
        module procedure cx_eq_int
        module procedure cx_eq_cx
        module procedure cx_eq_cdp
        module procedure cx_eq_csp
        module procedure cdp_eq_mf
        module procedure cdp_eq_cx
        module procedure csp_eq_mf
        module procedure csp_eq_cx
    end interface
    public :: operator(==)

    interface operator(/=)
        module procedure mf_ne_mf
        module procedure mf_ne_dp
        module procedure mf_ne_sp
        module procedure mf_ne_int
        module procedure mf_ne_cx
        module procedure mf_ne_cdp
        module procedure mf_ne_csp
        module procedure dp_ne_mf
        module procedure dp_ne_cx
        module procedure sp_ne_mf
        module procedure sp_ne_cx
        module procedure int_ne_mf
        module procedure int_ne_cx
        module procedure cx_ne_mf
        module procedure cx_ne_dp
        module procedure cx_ne_sp
        module procedure cx_ne_int
        module procedure cx_ne_cx
        module procedure cx_ne_cdp
        module procedure cx_ne_csp
        module procedure cdp_ne_mf
        module procedure cdp_ne_cx
        module procedure csp_ne_mf
        module procedure csp_ne_cx
    end interface
    public :: operator(/=)

    ! ================================================================
    ! Ordered comparison (<, >, <=, >=) — real types only
    ! ================================================================

    interface operator(<)
        module procedure mf_lt_mf
        module procedure mf_lt_dp
        module procedure mf_lt_sp
        module procedure mf_lt_int
        module procedure dp_lt_mf
        module procedure sp_lt_mf
        module procedure int_lt_mf
    end interface
    public :: operator(<)

    interface operator(>)
        module procedure mf_gt_mf
        module procedure mf_gt_dp
        module procedure mf_gt_sp
        module procedure mf_gt_int
        module procedure dp_gt_mf
        module procedure sp_gt_mf
        module procedure int_gt_mf
    end interface
    public :: operator(>)

    interface operator(<=)
        module procedure mf_le_mf
        module procedure mf_le_dp
        module procedure mf_le_sp
        module procedure mf_le_int
        module procedure dp_le_mf
        module procedure sp_le_mf
        module procedure int_le_mf
    end interface
    public :: operator(<=)

    interface operator(>=)
        module procedure mf_ge_mf
        module procedure mf_ge_dp
        module procedure mf_ge_sp
        module procedure mf_ge_int
        module procedure dp_ge_mf
        module procedure sp_ge_mf
        module procedure int_ge_mf
    end interface
    public :: operator(>=)

    ! ================================================================
    ! Unary math generics
    ! ================================================================

    interface abs
        module procedure mf_abs
        module procedure cx_abs
    end interface
    public :: abs

    interface acos
        module procedure mf_acos
        module procedure cx_acos
    end interface
    public :: acos

    interface aimag
        module procedure cx_aimag
    end interface
    public :: aimag

    interface asin
        module procedure mf_asin
        module procedure cx_asin
    end interface
    public :: asin

    interface atan
        module procedure mf_atan
        module procedure cx_atan
    end interface
    public :: atan

    interface conjg
        module procedure cx_conjg
    end interface
    public :: conjg

    interface cos
        module procedure mf_cos
        module procedure cx_cos
    end interface
    public :: cos

    interface exp
        module procedure mf_exp
        module procedure cx_exp
    end interface
    public :: exp

    interface log
        module procedure mf_log
        module procedure cx_log
    end interface
    public :: log

    interface log10
        module procedure mf_log10
    end interface
    public :: log10

    interface sin
        module procedure mf_sin
        module procedure cx_sin
    end interface
    public :: sin

    interface sqrt
        module procedure mf_sqrt
        module procedure cx_sqrt
    end interface
    public :: sqrt

    interface tan
        module procedure mf_tan
        module procedure cx_tan
    end interface
    public :: tan

    ! ================================================================
    ! Binary real functions (sign, mod, atan2, min, max)
    ! ================================================================

    interface sign
        module procedure mf_sign_mf
        module procedure mf_sign_dp
        module procedure mf_sign_sp
        module procedure mf_sign_int
        module procedure dp_sign_mf
        module procedure sp_sign_mf
        module procedure int_sign_mf
    end interface
    public :: sign

    interface mod
        module procedure mf_mod_mf
        module procedure mf_mod_dp
        module procedure mf_mod_sp
        module procedure mf_mod_int
        module procedure dp_mod_mf
        module procedure sp_mod_mf
        module procedure int_mod_mf
    end interface
    public :: mod

    interface atan2
        module procedure mf_atan2_mf
        module procedure mf_atan2_dp
        module procedure mf_atan2_sp
        module procedure mf_atan2_int
        module procedure dp_atan2_mf
        module procedure sp_atan2_mf
        module procedure int_atan2_mf
    end interface
    public :: atan2

    interface min
        module procedure mf_min_mf
        module procedure mf_min_dp
        module procedure mf_min_sp
        module procedure mf_min_int
        module procedure dp_min_mf
        module procedure sp_min_mf
        module procedure int_min_mf
        module procedure mf_min3
    end interface
    public :: min

    interface max
        module procedure mf_max_mf
        module procedure mf_max_dp
        module procedure mf_max_sp
        module procedure mf_max_int
        module procedure dp_max_mf
        module procedure sp_max_mf
        module procedure int_max_mf
        module procedure mf_max3
    end interface
    public :: max

    ! ================================================================
    ! Array reductions (maxval, minval, maxloc, minloc)
    ! ================================================================

    interface maxval
        module procedure mf_maxval_1d
        module procedure mf_maxval_1d_dim
        module procedure mf_maxval_1d_mask
        module procedure mf_maxval_1d_dim_mask
        module procedure mf_maxval_2d
        module procedure mf_maxval_2d_dim
        module procedure mf_maxval_2d_mask
        module procedure mf_maxval_2d_dim_mask
    end interface
    public :: maxval

    interface minval
        module procedure mf_minval_1d
        module procedure mf_minval_1d_dim
        module procedure mf_minval_1d_mask
        module procedure mf_minval_1d_dim_mask
        module procedure mf_minval_2d
        module procedure mf_minval_2d_dim
        module procedure mf_minval_2d_mask
        module procedure mf_minval_2d_dim_mask
    end interface
    public :: minval

    interface maxloc
        module procedure mf_maxloc_1d
        module procedure mf_maxloc_1d_dim
        module procedure mf_maxloc_1d_mask
        module procedure mf_maxloc_1d_dim_mask
        module procedure mf_maxloc_2d
        module procedure mf_maxloc_2d_dim
        module procedure mf_maxloc_2d_mask
        module procedure mf_maxloc_2d_dim_mask
    end interface
    public :: maxloc

    interface minloc
        module procedure mf_minloc_1d
        module procedure mf_minloc_1d_dim
        module procedure mf_minloc_1d_mask
        module procedure mf_minloc_1d_dim_mask
        module procedure mf_minloc_2d
        module procedure mf_minloc_2d_dim
        module procedure mf_minloc_2d_mask
        module procedure mf_minloc_2d_dim_mask
    end interface
    public :: minloc

    ! ================================================================
    ! Complex-specific: cmplx, real, mf_real
    ! ================================================================

    interface cmplx
        module procedure mf_cmplx_2
        module procedure mf_cmplx_1
    end interface
    public :: cmplx

    interface real
        module procedure cx_real
        module procedure mf_real_mf
    end interface
    public :: real

    interface mf_real
        module procedure mf_real_mf
        module procedure cx_real_mf
    end interface
    public :: mf_real

    ! ================================================================
    ! Conversions: dble, int, nint, mf_to_double
    ! ================================================================

    interface dble
        module procedure mf_dble
    end interface
    public :: dble

    interface int
        module procedure mf_int
    end interface
    public :: int

    interface nint
        module procedure mf_nint
    end interface
    public :: nint

    public :: mf_to_double

    ! ================================================================
    ! Assignment (=)
    ! ================================================================

    interface assignment(=)
        module procedure mf_assign_dp
        module procedure mf_assign_sp
        module procedure mf_assign_int
        module procedure mf_assign_cx
        module procedure mf_assign_cdp
        module procedure mf_assign_csp
        module procedure cx_assign_mf
        module procedure cx_assign_dp
        module procedure cx_assign_sp
        module procedure cx_assign_int
        module procedure cx_assign_cdp
        module procedure cx_assign_csp
    end interface
    public :: assignment(=)

    ! ================================================================
    ! Defined I/O
    ! ================================================================

    interface write(formatted)
        module procedure write_mf_formatted
        module procedure write_cx_formatted
    end interface
    public :: write(formatted)

    interface read(formatted)
        module procedure read_mf_formatted
        module procedure read_cx_formatted
    end interface
    public :: read(formatted)

    ! ================================================================
    ! Named constants
    ! ================================================================

    type(float64x2), parameter, public :: MF_ZERO  = float64x2([0.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_ONE   = float64x2([1.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_TWO   = float64x2([2.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_HALF  = float64x2([0.5_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_EIGHT = float64x2([8.0_dp, 0.0_dp])

    type(float64x2), parameter, public :: MF_SAFMIN = float64x2([tiny(0.0_dp),  0.0_dp])
    type(float64x2), parameter, public :: MF_SAFMAX = float64x2([huge(0.0_dp),  0.0_dp])
    type(float64x2), parameter, public :: MF_TSML   = float64x2([1.0e-100_dp,   0.0_dp])
    type(float64x2), parameter, public :: MF_TBIG   = float64x2([1.0e+100_dp,   0.0_dp])
    type(float64x2), parameter, public :: MF_SSML   = float64x2([1.0e-50_dp,    0.0_dp])
    type(float64x2), parameter, public :: MF_SBIG   = float64x2([1.0e+50_dp,    0.0_dp])
    type(float64x2), parameter, public :: MF_RTMIN  = float64x2([1.0e-150_dp,   0.0_dp])
    type(float64x2), parameter, public :: MF_RTMAX  = float64x2([1.0e+150_dp,   0.0_dp])

contains

    ! ================================================================
    ! Internal conversion helpers
    ! ================================================================

    pure function to_qp(mf) result(res)
        type(float64x2), intent(in) :: mf
        real(qp) :: res
        res = real(mf%limbs(1), qp) + real(mf%limbs(2), qp)
    end function

    pure function from_qp(v) result(res)
        real(qp), intent(in) :: v
        type(float64x2) :: res
        res%limbs(1) = real(v, dp)
        res%limbs(2) = real(v - real(res%limbs(1), qp), dp)
    end function

    pure function to_cqp(cx) result(res)
        type(complex128x2), intent(in) :: cx
        complex(qp) :: res
        res = cmplx(to_qp(cx%re), to_qp(cx%im), kind=qp)
    end function

    pure function from_cqp(v) result(res)
        complex(qp), intent(in) :: v
        type(complex128x2) :: res
        res%re = from_qp(real(v, qp))
        res%im = from_qp(aimag(v))
    end function

    ! ================================================================
    ! Constructors
    ! ================================================================

    pure function mf_from_dp(v) result(res)
        real(dp), intent(in) :: v
        type(float64x2) :: res
        res = from_qp(real(v, qp))
    end function

    pure function mf_from_sp(v) result(res)
        real(sp), intent(in) :: v
        type(float64x2) :: res
        res = from_qp(real(v, qp))
    end function

    pure function mf_from_int(v) result(res)
        integer, intent(in) :: v
        type(float64x2) :: res
        res = from_qp(real(v, qp))
    end function

    function mf_from_char(v) result(res)
        character(len=*), intent(in) :: v
        type(float64x2) :: res
        real(qp) :: val_qp
        read(v, *) val_qp
        res = from_qp(val_qp)
    end function

    pure function mf_from_cx(cx) result(res)
        type(complex128x2), intent(in) :: cx
        type(float64x2) :: res
        res = cx%re
    end function

    pure function mf_from_cdp(x) result(res)
        complex(dp), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(real(x, qp))
    end function

    pure function mf_from_csp(x) result(res)
        complex(sp), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(real(x, qp))
    end function

    pure function cx_from_mf_2(r, i) result(res)
        type(float64x2), intent(in) :: r, i
        type(complex128x2) :: res
        res%re = r
        res%im = i
    end function

    pure function cx_from_dp_2(r, i) result(res)
        real(dp), intent(in) :: r, i
        type(complex128x2) :: res
        res%re = from_qp(real(r, qp))
        res%im = from_qp(real(i, qp))
    end function

    pure function cx_from_mf_1(r) result(res)
        type(float64x2), intent(in) :: r
        type(complex128x2) :: res
        res%re = r
        res%im = from_qp(0.0_qp)
    end function

    pure function cx_from_dp_1(r) result(res)
        real(dp), intent(in) :: r
        type(complex128x2) :: res
        res%re = from_qp(real(r, qp))
        res%im = from_qp(0.0_qp)
    end function

    ! ================================================================
    ! Assignment
    ! ================================================================

    pure subroutine mf_assign_dp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        real(dp), intent(in) :: rhs
        lhs = from_qp(real(rhs, qp))
    end subroutine

    pure subroutine mf_assign_sp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        real(sp), intent(in) :: rhs
        lhs = from_qp(real(rhs, qp))
    end subroutine

    pure subroutine mf_assign_int(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        integer, intent(in) :: rhs
        lhs = from_qp(real(rhs, qp))
    end subroutine

    pure subroutine mf_assign_cx(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        type(complex128x2), intent(in) :: rhs
        lhs = from_qp(real(to_cqp(rhs), qp))
    end subroutine

    pure subroutine mf_assign_cdp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        complex(dp), intent(in) :: rhs
        lhs = from_qp(real(cmplx(rhs, kind=qp), qp))
    end subroutine

    pure subroutine mf_assign_csp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        complex(sp), intent(in) :: rhs
        lhs = from_qp(real(cmplx(rhs, kind=qp), qp))
    end subroutine

    pure subroutine cx_assign_mf(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        type(float64x2), intent(in) :: rhs
        lhs = from_cqp(cmplx(to_qp(rhs), 0.0_qp, qp))
    end subroutine

    pure subroutine cx_assign_dp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        real(dp), intent(in) :: rhs
        lhs = from_cqp(cmplx(real(rhs, qp), 0.0_qp, qp))
    end subroutine

    pure subroutine cx_assign_sp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        real(sp), intent(in) :: rhs
        lhs = from_cqp(cmplx(real(rhs, qp), 0.0_qp, qp))
    end subroutine

    pure subroutine cx_assign_int(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        integer, intent(in) :: rhs
        lhs = from_cqp(cmplx(real(rhs, qp), 0.0_qp, qp))
    end subroutine

    pure subroutine cx_assign_cdp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        complex(dp), intent(in) :: rhs
        lhs = from_cqp(cmplx(rhs, kind=qp))
    end subroutine

    pure subroutine cx_assign_csp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        complex(sp), intent(in) :: rhs
        lhs = from_cqp(cmplx(rhs, kind=qp))
    end subroutine

    ! ================================================================
    ! Conversions
    ! ================================================================

    pure function mf_dble(x) result(res)
        type(float64x2), intent(in) :: x
        real(dp) :: res
        res = real(to_qp(x), dp)
    end function

    pure function mf_int(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = int(to_qp(x))
    end function

    pure function mf_nint(x) result(res)
        type(float64x2), intent(in) :: x
        integer :: res
        res = nint(to_qp(x))
    end function

    pure function mf_to_double(x) result(res)
        type(float64x2), intent(in) :: x
        real(dp) :: res
        res = real(to_qp(x), dp)
    end function

    ! ================================================================
    ! Unary operators
    ! ================================================================

    pure function mf_neg(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(-to_qp(x))
    end function

    pure function cx_neg(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(-to_cqp(x))
    end function

    ! ================================================================
    ! Unary math: mf → mf
    ! ================================================================

    pure function mf_abs(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(abs(to_qp(x)))
    end function

    pure function mf_sqrt(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(sqrt(to_qp(x)))
    end function

    pure function mf_sin(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(sin(to_qp(x)))
    end function

    pure function mf_cos(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(cos(to_qp(x)))
    end function

    pure function mf_tan(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(tan(to_qp(x)))
    end function

    pure function mf_exp(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(exp(to_qp(x)))
    end function

    pure function mf_log(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(log(to_qp(x)))
    end function

    pure function mf_log10(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(log10(to_qp(x)))
    end function

    pure function mf_atan(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(atan(to_qp(x)))
    end function

    pure function mf_asin(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(asin(to_qp(x)))
    end function

    pure function mf_acos(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(acos(to_qp(x)))
    end function

    ! ================================================================
    ! Unary math: cx → cx
    ! ================================================================

    pure function cx_sqrt(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(sqrt(to_cqp(x)))
    end function

    pure function cx_sin(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(sin(to_cqp(x)))
    end function

    pure function cx_cos(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(cos(to_cqp(x)))
    end function

    pure function cx_tan(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(tan(to_cqp(x)))
    end function

    pure function cx_exp(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(exp(to_cqp(x)))
    end function

    pure function cx_log(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(log(to_cqp(x)))
    end function

    pure function cx_atan(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(atan(to_cqp(x)))
    end function

    pure function cx_asin(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(asin(to_cqp(x)))
    end function

    pure function cx_acos(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(acos(to_cqp(x)))
    end function

    pure function cx_conjg(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(conjg(to_cqp(x)))
    end function

    ! ================================================================
    ! Unary math: cx → mf (abs, aimag)
    ! ================================================================

    pure function cx_abs(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(abs(to_cqp(x)))
    end function

    pure function cx_aimag(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(aimag(to_cqp(x)))
    end function

    ! ================================================================
    ! Complex-specific: cmplx, real, mf_real
    ! ================================================================

    pure function mf_cmplx_2(re, im) result(res)
        type(float64x2), intent(in) :: re, im
        type(complex128x2) :: res
        res%re = re
        res%im = im
    end function

    pure function mf_cmplx_1(re) result(res)
        type(float64x2), intent(in) :: re
        type(complex128x2) :: res
        res%re = re
        res%im = from_qp(0.0_qp)
    end function

    pure function cx_real(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(real(to_cqp(x), qp))
    end function

    pure function mf_real_mf(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res = x
    end function

    pure function cx_real_mf(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = x%re
    end function

    ! ================================================================
    ! Binary real functions (sign, mod, atan2, min, max)
    ! ================================================================

    pure function mf_sign_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(to_qp(a), to_qp(b)))
    end function

    pure function mf_sign_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(to_qp(a), real(b, qp)))
    end function

    pure function mf_sign_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(to_qp(a), real(b, qp)))
    end function

    pure function mf_sign_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(to_qp(a), real(b, qp)))
    end function

    pure function dp_sign_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(real(a, qp), to_qp(b)))
    end function

    pure function sp_sign_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(real(a, qp), to_qp(b)))
    end function

    pure function int_sign_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(sign(real(a, qp), to_qp(b)))
    end function

    pure function mf_mod_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(to_qp(a), to_qp(b)))
    end function

    pure function mf_mod_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(to_qp(a), real(b, qp)))
    end function

    pure function mf_mod_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(to_qp(a), real(b, qp)))
    end function

    pure function mf_mod_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(to_qp(a), real(b, qp)))
    end function

    pure function dp_mod_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(real(a, qp), to_qp(b)))
    end function

    pure function sp_mod_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(real(a, qp), to_qp(b)))
    end function

    pure function int_mod_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(mod(real(a, qp), to_qp(b)))
    end function

    pure function mf_atan2_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(a), to_qp(b)))
    end function

    pure function mf_atan2_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(a), real(b, qp)))
    end function

    pure function mf_atan2_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(a), real(b, qp)))
    end function

    pure function mf_atan2_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(a), real(b, qp)))
    end function

    pure function dp_atan2_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(real(a, qp), to_qp(b)))
    end function

    pure function sp_atan2_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(real(a, qp), to_qp(b)))
    end function

    pure function int_atan2_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(atan2(real(a, qp), to_qp(b)))
    end function

    pure function mf_min_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), to_qp(b)))
    end function

    pure function mf_min_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), real(b, qp)))
    end function

    pure function mf_min_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), real(b, qp)))
    end function

    pure function mf_min_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), real(b, qp)))
    end function

    pure function dp_min_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(real(a, qp), to_qp(b)))
    end function

    pure function sp_min_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(real(a, qp), to_qp(b)))
    end function

    pure function int_min_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(min(real(a, qp), to_qp(b)))
    end function

    pure function mf_max_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), to_qp(b)))
    end function

    pure function mf_max_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), real(b, qp)))
    end function

    pure function mf_max_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), real(b, qp)))
    end function

    pure function mf_max_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), real(b, qp)))
    end function

    pure function dp_max_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(real(a, qp), to_qp(b)))
    end function

    pure function sp_max_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(real(a, qp), to_qp(b)))
    end function

    pure function int_max_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(max(real(a, qp), to_qp(b)))
    end function

    ! min/max with 3 arguments
    pure function mf_min3(a, b, c) result(res)
        type(float64x2), intent(in) :: a, b, c
        type(float64x2) :: res
        res = from_qp(min(to_qp(a), to_qp(b), to_qp(c)))
    end function

    pure function mf_max3(a, b, c) result(res)
        type(float64x2), intent(in) :: a, b, c
        type(float64x2) :: res
        res = from_qp(max(to_qp(a), to_qp(b), to_qp(c)))
    end function

    ! ================================================================
    ! Array reductions (maxval, minval, maxloc, minloc)
    ! ================================================================

    ! maxval — rank 1
    pure function mf_maxval_1d(array) result(res)
        type(float64x2), intent(in) :: array(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = from_qp(maxval(tmp))
    end function

    pure function mf_maxval_1d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = from_qp(maxval(tmp, dim=dim))
    end function

    pure function mf_maxval_1d_mask(array, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in) :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = from_qp(maxval(tmp, mask=mask))
    end function

    pure function mf_maxval_1d_dim_mask(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in) :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = from_qp(maxval(tmp, dim=dim, mask=mask))
    end function

    ! maxval — rank 2
    pure function mf_maxval_2d(array) result(res)
        type(float64x2), intent(in) :: array(:,:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = from_qp(maxval(tmp))
    end function

    pure function mf_maxval_2d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer, intent(in) :: dim
        type(float64x2) :: res(merge(size(array,2), size(array,1), dim==1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2), size(array,1), dim==1))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        tmp_res = maxval(tmp, dim=dim)
        do i = 1, size(tmp_res)
            res(i) = from_qp(tmp_res(i))
        end do
    end function

    pure function mf_maxval_2d_mask(array, mask) result(res)
        type(float64x2), intent(in) :: array(:,:)
        logical, intent(in) :: mask(:,:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = from_qp(maxval(tmp, mask=mask))
    end function

    pure function mf_maxval_2d_dim_mask(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer, intent(in) :: dim
        logical, intent(in) :: mask(:,:)
        type(float64x2) :: res(merge(size(array,2), size(array,1), dim==1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2), size(array,1), dim==1))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        tmp_res = maxval(tmp, dim=dim, mask=mask)
        do i = 1, size(tmp_res)
            res(i) = from_qp(tmp_res(i))
        end do
    end function

    ! minval — rank 1
    pure function mf_minval_1d(array) result(res)
        type(float64x2), intent(in) :: array(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = from_qp(minval(tmp))
    end function

    pure function mf_minval_1d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = from_qp(minval(tmp, dim=dim))
    end function

    pure function mf_minval_1d_mask(array, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in) :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = from_qp(minval(tmp, mask=mask))
    end function

    pure function mf_minval_1d_dim_mask(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in) :: mask(:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = from_qp(minval(tmp, dim=dim, mask=mask))
    end function

    ! minval — rank 2
    pure function mf_minval_2d(array) result(res)
        type(float64x2), intent(in) :: array(:,:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = from_qp(minval(tmp))
    end function

    pure function mf_minval_2d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer, intent(in) :: dim
        type(float64x2) :: res(merge(size(array,2), size(array,1), dim==1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2), size(array,1), dim==1))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        tmp_res = minval(tmp, dim=dim)
        do i = 1, size(tmp_res)
            res(i) = from_qp(tmp_res(i))
        end do
    end function

    pure function mf_minval_2d_mask(array, mask) result(res)
        type(float64x2), intent(in) :: array(:,:)
        logical, intent(in) :: mask(:,:)
        type(float64x2) :: res
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = from_qp(minval(tmp, mask=mask))
    end function

    pure function mf_minval_2d_dim_mask(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer, intent(in) :: dim
        logical, intent(in) :: mask(:,:)
        type(float64x2) :: res(merge(size(array,2), size(array,1), dim==1))
        real(qp) :: tmp(size(array,1), size(array,2))
        real(qp) :: tmp_res(merge(size(array,2), size(array,1), dim==1))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        tmp_res = minval(tmp, dim=dim, mask=mask)
        do i = 1, size(tmp_res)
            res(i) = from_qp(tmp_res(i))
        end do
    end function

    ! maxloc — rank 1
    pure function mf_maxloc_1d(array) result(res)
        type(float64x2), intent(in) :: array(:)
        integer :: res(1)
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = maxloc(tmp)
    end function

    pure function mf_maxloc_1d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        integer :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = maxloc(tmp, dim=dim)
    end function

    pure function mf_maxloc_1d_mask(array, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in) :: mask(:)
        integer :: res(1)
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = maxloc(tmp, mask=mask)
    end function

    pure function mf_maxloc_1d_dim_mask(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in) :: mask(:)
        integer :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = maxloc(tmp, dim=dim, mask=mask)
    end function

    ! maxloc — rank 2
    pure function mf_maxloc_2d(array) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer :: res(2)
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = maxloc(tmp)
    end function

    pure function mf_maxloc_2d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer, intent(in) :: dim
        integer :: res(merge(size(array,2), size(array,1), dim==1))
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = maxloc(tmp, dim=dim)
    end function

    pure function mf_maxloc_2d_mask(array, mask) result(res)
        type(float64x2), intent(in) :: array(:,:)
        logical, intent(in) :: mask(:,:)
        integer :: res(2)
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = maxloc(tmp, mask=mask)
    end function

    pure function mf_maxloc_2d_dim_mask(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer, intent(in) :: dim
        logical, intent(in) :: mask(:,:)
        integer :: res(merge(size(array,2), size(array,1), dim==1))
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = maxloc(tmp, dim=dim, mask=mask)
    end function

    ! minloc — rank 1
    pure function mf_minloc_1d(array) result(res)
        type(float64x2), intent(in) :: array(:)
        integer :: res(1)
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = minloc(tmp)
    end function

    pure function mf_minloc_1d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        integer :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = minloc(tmp, dim=dim)
    end function

    pure function mf_minloc_1d_mask(array, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        logical, intent(in) :: mask(:)
        integer :: res(1)
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = minloc(tmp, mask=mask)
    end function

    pure function mf_minloc_1d_dim_mask(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:)
        integer, intent(in) :: dim
        logical, intent(in) :: mask(:)
        integer :: res
        real(qp) :: tmp(size(array))
        integer :: i
        do i = 1, size(array)
            tmp(i) = to_qp(array(i))
        end do
        res = minloc(tmp, dim=dim, mask=mask)
    end function

    ! minloc — rank 2
    pure function mf_minloc_2d(array) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer :: res(2)
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = minloc(tmp)
    end function

    pure function mf_minloc_2d_dim(array, dim) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer, intent(in) :: dim
        integer :: res(merge(size(array,2), size(array,1), dim==1))
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = minloc(tmp, dim=dim)
    end function

    pure function mf_minloc_2d_mask(array, mask) result(res)
        type(float64x2), intent(in) :: array(:,:)
        logical, intent(in) :: mask(:,:)
        integer :: res(2)
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = minloc(tmp, mask=mask)
    end function

    pure function mf_minloc_2d_dim_mask(array, dim, mask) result(res)
        type(float64x2), intent(in) :: array(:,:)
        integer, intent(in) :: dim
        logical, intent(in) :: mask(:,:)
        integer :: res(merge(size(array,2), size(array,1), dim==1))
        real(qp) :: tmp(size(array,1), size(array,2))
        integer :: i, j
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                tmp(i,j) = to_qp(array(i,j))
            end do
        end do
        res = minloc(tmp, dim=dim, mask=mask)
    end function

    ! ================================================================
    ! Binary arithmetic (+, -, *, /)
    ! ================================================================

    pure function mf_add_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + to_qp(b))
    end function

    pure function mf_add_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + real(b, qp))
    end function

    pure function mf_add_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + real(b, qp))
    end function

    pure function mf_add_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + real(b, qp))
    end function

    pure function mf_add_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) + to_cqp(b))
    end function

    pure function mf_add_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) + cmplx(b, kind=qp))
    end function

    pure function mf_add_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) + cmplx(b, kind=qp))
    end function

    pure function dp_add_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function

    pure function dp_add_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) + to_cqp(b))
    end function

    pure function sp_add_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function

    pure function sp_add_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) + to_cqp(b))
    end function

    pure function int_add_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function

    pure function int_add_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) + to_cqp(b))
    end function

    pure function cx_add_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + to_qp(b))
    end function

    pure function cx_add_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + real(b, qp))
    end function

    pure function cx_add_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + real(b, qp))
    end function

    pure function cx_add_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + real(b, qp))
    end function

    pure function cx_add_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + to_cqp(b))
    end function

    pure function cx_add_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + cmplx(b, kind=qp))
    end function

    pure function cx_add_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + cmplx(b, kind=qp))
    end function

    pure function cdp_add_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) + to_qp(b))
    end function

    pure function cdp_add_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) + to_cqp(b))
    end function

    pure function csp_add_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) + to_qp(b))
    end function

    pure function csp_add_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) + to_cqp(b))
    end function

    pure function mf_sub_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - to_qp(b))
    end function

    pure function mf_sub_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - real(b, qp))
    end function

    pure function mf_sub_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - real(b, qp))
    end function

    pure function mf_sub_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - real(b, qp))
    end function

    pure function mf_sub_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) - to_cqp(b))
    end function

    pure function mf_sub_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) - cmplx(b, kind=qp))
    end function

    pure function mf_sub_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) - cmplx(b, kind=qp))
    end function

    pure function dp_sub_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function

    pure function dp_sub_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) - to_cqp(b))
    end function

    pure function sp_sub_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function

    pure function sp_sub_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) - to_cqp(b))
    end function

    pure function int_sub_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function

    pure function int_sub_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) - to_cqp(b))
    end function

    pure function cx_sub_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - to_qp(b))
    end function

    pure function cx_sub_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - real(b, qp))
    end function

    pure function cx_sub_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - real(b, qp))
    end function

    pure function cx_sub_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - real(b, qp))
    end function

    pure function cx_sub_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - to_cqp(b))
    end function

    pure function cx_sub_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - cmplx(b, kind=qp))
    end function

    pure function cx_sub_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - cmplx(b, kind=qp))
    end function

    pure function cdp_sub_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) - to_qp(b))
    end function

    pure function cdp_sub_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) - to_cqp(b))
    end function

    pure function csp_sub_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) - to_qp(b))
    end function

    pure function csp_sub_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) - to_cqp(b))
    end function

    pure function mf_mul_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * to_qp(b))
    end function

    pure function mf_mul_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * real(b, qp))
    end function

    pure function mf_mul_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * real(b, qp))
    end function

    pure function mf_mul_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * real(b, qp))
    end function

    pure function mf_mul_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) * to_cqp(b))
    end function

    pure function mf_mul_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) * cmplx(b, kind=qp))
    end function

    pure function mf_mul_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) * cmplx(b, kind=qp))
    end function

    pure function dp_mul_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function

    pure function dp_mul_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) * to_cqp(b))
    end function

    pure function sp_mul_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function

    pure function sp_mul_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) * to_cqp(b))
    end function

    pure function int_mul_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function

    pure function int_mul_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) * to_cqp(b))
    end function

    pure function cx_mul_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * to_qp(b))
    end function

    pure function cx_mul_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * real(b, qp))
    end function

    pure function cx_mul_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * real(b, qp))
    end function

    pure function cx_mul_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * real(b, qp))
    end function

    pure function cx_mul_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * to_cqp(b))
    end function

    pure function cx_mul_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * cmplx(b, kind=qp))
    end function

    pure function cx_mul_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * cmplx(b, kind=qp))
    end function

    pure function cdp_mul_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) * to_qp(b))
    end function

    pure function cdp_mul_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) * to_cqp(b))
    end function

    pure function csp_mul_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) * to_qp(b))
    end function

    pure function csp_mul_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) * to_cqp(b))
    end function

    pure function mf_div_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / to_qp(b))
    end function

    pure function mf_div_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / real(b, qp))
    end function

    pure function mf_div_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / real(b, qp))
    end function

    pure function mf_div_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / real(b, qp))
    end function

    pure function mf_div_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) / to_cqp(b))
    end function

    pure function mf_div_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) / cmplx(b, kind=qp))
    end function

    pure function mf_div_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) / cmplx(b, kind=qp))
    end function

    pure function dp_div_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function

    pure function dp_div_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) / to_cqp(b))
    end function

    pure function sp_div_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function

    pure function sp_div_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) / to_cqp(b))
    end function

    pure function int_div_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function

    pure function int_div_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) / to_cqp(b))
    end function

    pure function cx_div_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / to_qp(b))
    end function

    pure function cx_div_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / real(b, qp))
    end function

    pure function cx_div_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / real(b, qp))
    end function

    pure function cx_div_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / real(b, qp))
    end function

    pure function cx_div_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / to_cqp(b))
    end function

    pure function cx_div_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / cmplx(b, kind=qp))
    end function

    pure function cx_div_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / cmplx(b, kind=qp))
    end function

    pure function cdp_div_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) / to_qp(b))
    end function

    pure function cdp_div_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) / to_cqp(b))
    end function

    pure function csp_div_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) / to_qp(b))
    end function

    pure function csp_div_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) / to_cqp(b))
    end function

    ! ================================================================
    ! Power operator (**)  — preserves integer exponent semantics
    ! ================================================================

    pure function mf_pow_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** to_qp(b))
    end function

    pure function mf_pow_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** real(b, qp))
    end function

    pure function mf_pow_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** real(b, qp))
    end function

    pure function mf_pow_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** b)
    end function

    pure function mf_pow_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) ** to_cqp(b))
    end function

    pure function mf_pow_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) ** cmplx(b, kind=qp))
    end function

    pure function mf_pow_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) ** cmplx(b, kind=qp))
    end function

    pure function dp_pow_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) ** to_qp(b))
    end function

    pure function dp_pow_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) ** to_cqp(b))
    end function

    pure function sp_pow_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) ** to_qp(b))
    end function

    pure function sp_pow_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) ** to_cqp(b))
    end function

    pure function int_pow_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) ** to_qp(b))
    end function

    pure function int_pow_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(real(a, qp) ** to_cqp(b))
    end function

    pure function cx_pow_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** to_qp(b))
    end function

    pure function cx_pow_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** real(b, qp))
    end function

    pure function cx_pow_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** real(b, qp))
    end function

    pure function cx_pow_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** b)
    end function

    pure function cx_pow_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** to_cqp(b))
    end function

    pure function cx_pow_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** cmplx(b, kind=qp))
    end function

    pure function cx_pow_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) ** cmplx(b, kind=qp))
    end function

    pure function cdp_pow_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) ** to_qp(b))
    end function

    pure function cdp_pow_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) ** to_cqp(b))
    end function

    pure function csp_pow_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) ** to_qp(b))
    end function

    pure function csp_pow_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(cmplx(a, kind=qp) ** to_cqp(b))
    end function

    ! ================================================================
    ! Equality comparison (==, /=)
    ! ================================================================

    pure function mf_eq_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) == to_qp(b)
    end function

    pure function mf_eq_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function

    pure function mf_eq_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function

    pure function mf_eq_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function

    pure function mf_eq_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_qp(a) == to_cqp(b)
    end function

    pure function mf_eq_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) == cmplx(b, kind=qp)
    end function

    pure function mf_eq_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) == cmplx(b, kind=qp)
    end function

    pure function dp_eq_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function

    pure function dp_eq_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_cqp(b)
    end function

    pure function sp_eq_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function

    pure function sp_eq_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_cqp(b)
    end function

    pure function int_eq_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function

    pure function int_eq_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_cqp(b)
    end function

    pure function cx_eq_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) == to_qp(b)
    end function

    pure function cx_eq_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_cqp(a) == real(b, qp)
    end function

    pure function cx_eq_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_cqp(a) == real(b, qp)
    end function

    pure function cx_eq_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_cqp(a) == real(b, qp)
    end function

    pure function cx_eq_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) == to_cqp(b)
    end function

    pure function cx_eq_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        logical :: res
        res = to_cqp(a) == cmplx(b, kind=qp)
    end function

    pure function cx_eq_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        logical :: res
        res = to_cqp(a) == cmplx(b, kind=qp)
    end function

    pure function cdp_eq_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) == to_qp(b)
    end function

    pure function cdp_eq_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) == to_cqp(b)
    end function

    pure function csp_eq_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) == to_qp(b)
    end function

    pure function csp_eq_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) == to_cqp(b)
    end function

    pure function mf_ne_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) /= to_qp(b)
    end function

    pure function mf_ne_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function

    pure function mf_ne_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function

    pure function mf_ne_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function

    pure function mf_ne_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_qp(a) /= to_cqp(b)
    end function

    pure function mf_ne_cdp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= cmplx(b, kind=qp)
    end function

    pure function mf_ne_csp(a, b) result(res)
        type(float64x2), intent(in) :: a
        complex(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= cmplx(b, kind=qp)
    end function

    pure function dp_ne_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function

    pure function dp_ne_cx(a, b) result(res)
        real(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_cqp(b)
    end function

    pure function sp_ne_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function

    pure function sp_ne_cx(a, b) result(res)
        real(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_cqp(b)
    end function

    pure function int_ne_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function

    pure function int_ne_cx(a, b) result(res)
        integer, intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_cqp(b)
    end function

    pure function cx_ne_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= to_qp(b)
    end function

    pure function cx_ne_dp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= real(b, qp)
    end function

    pure function cx_ne_sp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= real(b, qp)
    end function

    pure function cx_ne_int(a, b) result(res)
        type(complex128x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_cqp(a) /= real(b, qp)
    end function

    pure function cx_ne_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= to_cqp(b)
    end function

    pure function cx_ne_cdp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(dp), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= cmplx(b, kind=qp)
    end function

    pure function cx_ne_csp(a, b) result(res)
        type(complex128x2), intent(in) :: a
        complex(sp), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= cmplx(b, kind=qp)
    end function

    pure function cdp_ne_mf(a, b) result(res)
        complex(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) /= to_qp(b)
    end function

    pure function cdp_ne_cx(a, b) result(res)
        complex(dp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) /= to_cqp(b)
    end function

    pure function csp_ne_mf(a, b) result(res)
        complex(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) /= to_qp(b)
    end function

    pure function csp_ne_cx(a, b) result(res)
        complex(sp), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = cmplx(a, kind=qp) /= to_cqp(b)
    end function

    ! ================================================================
    ! Ordered comparison (<, >, <=, >=)
    ! ================================================================

    pure function mf_lt_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) < to_qp(b)
    end function

    pure function mf_lt_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function

    pure function mf_lt_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function

    pure function mf_lt_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function

    pure function dp_lt_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function

    pure function sp_lt_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function

    pure function int_lt_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function

    pure function mf_gt_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) > to_qp(b)
    end function

    pure function mf_gt_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function

    pure function mf_gt_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function

    pure function mf_gt_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function

    pure function dp_gt_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function

    pure function sp_gt_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function

    pure function int_gt_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function

    pure function mf_le_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) <= to_qp(b)
    end function

    pure function mf_le_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function

    pure function mf_le_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function

    pure function mf_le_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function

    pure function dp_le_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function

    pure function sp_le_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function

    pure function int_le_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function

    pure function mf_ge_mf(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_qp(a) >= to_qp(b)
    end function

    pure function mf_ge_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function

    pure function mf_ge_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function

    pure function mf_ge_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function

    pure function dp_ge_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function

    pure function sp_ge_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function

    pure function int_ge_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function

    ! ================================================================
    ! Defined I/O
    ! ================================================================

    subroutine write_mf_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
        class(float64x2), intent(in) :: dtv
        integer, intent(in) :: unit
        character(*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(*), intent(inout) :: iomsg
        real(qp) :: val_qp
        val_qp = to_qp(dtv)
        if (iotype == 'LISTDIRECTED') then
            write(unit, *, iostat=iostat, iomsg=iomsg) val_qp
        else if (iotype(1:2) == 'DT' .and. len(iotype) > 2) then
            write(unit, fmt='(' // iotype(3:) // ')', iostat=iostat, iomsg=iomsg) val_qp
        else
            write(unit, *, iostat=iostat, iomsg=iomsg) val_qp
        end if
    end subroutine

    subroutine write_cx_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
        class(complex128x2), intent(in) :: dtv
        integer, intent(in) :: unit
        character(*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(*), intent(inout) :: iomsg
        complex(qp) :: val_cqp
        val_cqp = to_cqp(dtv)
        if (iotype == 'LISTDIRECTED') then
            write(unit, *, iostat=iostat, iomsg=iomsg) val_cqp
        else if (iotype(1:2) == 'DT' .and. len(iotype) > 2) then
            write(unit, fmt='(' // iotype(3:) // ')', iostat=iostat, iomsg=iomsg) val_cqp
        else
            write(unit, *, iostat=iostat, iomsg=iomsg) val_cqp
        end if
    end subroutine

    subroutine read_mf_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
        class(float64x2), intent(inout) :: dtv
        integer, intent(in) :: unit
        character(*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(*), intent(inout) :: iomsg
        real(qp) :: val_qp
        if (iotype == 'LISTDIRECTED') then
            read(unit, *, iostat=iostat, iomsg=iomsg) val_qp
        else if (iotype(1:2) == 'DT' .and. len(iotype) > 2) then
            read(unit, fmt='(' // iotype(3:) // ')', iostat=iostat, iomsg=iomsg) val_qp
        else
            read(unit, *, iostat=iostat, iomsg=iomsg) val_qp
        end if
        if (iostat == 0) then
            dtv%limbs(1) = real(val_qp, dp)
            dtv%limbs(2) = real(val_qp - real(dtv%limbs(1), qp), dp)
        end if
    end subroutine

    subroutine read_cx_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
        class(complex128x2), intent(inout) :: dtv
        integer, intent(in) :: unit
        character(*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(*), intent(inout) :: iomsg
        complex(qp) :: val_cqp
        if (iotype == 'LISTDIRECTED') then
            read(unit, *, iostat=iostat, iomsg=iomsg) val_cqp
        else if (iotype(1:2) == 'DT' .and. len(iotype) > 2) then
            read(unit, fmt='(' // iotype(3:) // ')', iostat=iostat, iomsg=iomsg) val_cqp
        else
            read(unit, *, iostat=iostat, iomsg=iomsg) val_cqp
        end if
        if (iostat == 0) then
            dtv%re%limbs(1) = real(real(val_cqp, qp), dp)
            dtv%re%limbs(2) = real(real(val_cqp, qp) - real(dtv%re%limbs(1), qp), dp)
            dtv%im%limbs(1) = real(aimag(val_cqp), dp)
            dtv%im%limbs(2) = real(aimag(val_cqp) - real(dtv%im%limbs(1), qp), dp)
        end if
    end subroutine

end module multifloats
