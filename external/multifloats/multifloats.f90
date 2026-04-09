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

    ! Constructors
    interface float64x2
        module procedure mf_from_dp
        module procedure mf_from_sp
        module procedure mf_from_int
        module procedure mf_from_char
        module procedure mf_from_cx
        module procedure mf_from_cdp
    end interface float64x2

    interface complex128x2
        module procedure cx_from_mf
        module procedure cx_from_dp
        module procedure cx_from_mf_1
        module procedure cx_from_dp_1
    end interface complex128x2

    ! Arithmetic Operators
    interface operator(+)
        module procedure mf_add_mf
        module procedure mf_add_dp
        module procedure dp_add_mf
        module procedure mf_add_sp
        module procedure sp_add_mf
        module procedure mf_add_int
        module procedure int_add_mf
        module procedure cx_add_cx
        module procedure cx_add_mf
        module procedure mf_add_cx
    end interface
    public :: operator(+)

    interface operator(-)
        module procedure mf_sub_mf
        module procedure mf_sub_dp
        module procedure dp_sub_mf
        module procedure mf_sub_sp
        module procedure sp_sub_mf
        module procedure mf_sub_int
        module procedure int_sub_mf
        module procedure mf_neg
        module procedure cx_sub_cx
        module procedure cx_sub_mf
        module procedure mf_sub_cx
        module procedure cx_neg
    end interface
    public :: operator(-)

    interface operator(*)
        module procedure mf_mul_mf
        module procedure mf_mul_dp
        module procedure dp_mul_mf
        module procedure mf_mul_sp
        module procedure sp_mul_mf
        module procedure mf_mul_int
        module procedure int_mul_mf
        module procedure cx_mul_cx
        module procedure cx_mul_mf
        module procedure mf_mul_cx
    end interface
    public :: operator(*)

    interface operator(/)
        module procedure mf_div_mf
        module procedure mf_div_dp
        module procedure dp_div_mf
        module procedure mf_div_sp
        module procedure sp_div_mf
        module procedure mf_div_int
        module procedure int_div_mf
        module procedure cx_div_cx
        module procedure cx_div_mf
        module procedure mf_div_cx
    end interface
    public :: operator(/)

    interface operator(**)
        module procedure mf_pow_int
        module procedure mf_pow_mf
    end interface
    public :: operator(**)

    ! Comparison Operators
    interface operator(==)
        module procedure mf_eq_mf
        module procedure mf_eq_dp
        module procedure dp_eq_mf
        module procedure mf_eq_sp
        module procedure sp_eq_mf
        module procedure mf_eq_int
        module procedure int_eq_mf
        module procedure cx_eq_cx
        module procedure cx_eq_mf
        module procedure mf_eq_cx
    end interface
    public :: operator(==)

    interface operator(/=)
        module procedure mf_ne_mf
        module procedure mf_ne_dp
        module procedure dp_ne_mf
        module procedure mf_ne_sp
        module procedure sp_ne_mf
        module procedure mf_ne_int
        module procedure int_ne_mf
        module procedure cx_ne_cx
        module procedure cx_ne_mf
        module procedure mf_ne_cx
    end interface
    public :: operator(/=)

    interface operator(<)
        module procedure mf_lt_mf
        module procedure mf_lt_dp
        module procedure dp_lt_mf
        module procedure mf_lt_sp
        module procedure sp_lt_mf
        module procedure mf_lt_int
        module procedure int_lt_mf
    end interface
    public :: operator(<)

    interface operator(>)
        module procedure mf_gt_mf
        module procedure mf_gt_dp
        module procedure dp_gt_mf
        module procedure mf_gt_sp
        module procedure sp_gt_mf
        module procedure mf_gt_int
        module procedure int_gt_mf
    end interface
    public :: operator(>)

    interface operator(<=)
        module procedure mf_le_mf
        module procedure mf_le_dp
        module procedure dp_le_mf
        module procedure mf_le_sp
        module procedure sp_le_mf
        module procedure mf_le_int
        module procedure int_le_mf
    end interface
    public :: operator(<=)

    interface operator(>=)
        module procedure mf_ge_mf
        module procedure mf_ge_dp
        module procedure dp_ge_mf
        module procedure mf_ge_sp
        module procedure sp_ge_mf
        module procedure mf_ge_int
        module procedure int_ge_mf
    end interface
    public :: operator(>=)

    ! Math Generics
    interface abs
        module procedure mf_abs
        module procedure cx_abs
    end interface
    public :: abs

    interface sqrt
        module procedure mf_sqrt
        module procedure cx_sqrt
    end interface
    public :: sqrt

    interface sin
        module procedure mf_sin
        module procedure cx_sin
    end interface
    public :: sin

    interface cos
        module procedure mf_cos
        module procedure cx_cos
    end interface
    public :: cos

    interface tan
        module procedure mf_tan
        module procedure cx_tan
    end interface
    public :: tan

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

    interface atan2
        module procedure mf_atan2
    end interface
    public :: atan2

    interface atan
        module procedure mf_atan
        module procedure cx_atan
    end interface
    public :: atan

    interface asin
        module procedure mf_asin
        module procedure cx_asin
    end interface
    public :: asin

    interface acos
        module procedure mf_acos
        module procedure cx_acos
    end interface
    public :: acos

    interface aimag
        module procedure cx_aimag
    end interface
    public :: aimag

    interface conjg
        module procedure cx_conjg
    end interface
    public :: conjg

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

    interface min
        module procedure mf_min
        module procedure mf_min3
    end interface
    public :: min

    interface max
        module procedure mf_max
        module procedure mf_max3
    end interface
    public :: max

    interface sign
        module procedure mf_sign
    end interface
    public :: sign

    interface mod
        module procedure mf_mod
    end interface
    public :: mod

    ! Conversion
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

    interface assignment(=)
        module procedure mf_assign_dp
        module procedure mf_assign_sp
        module procedure mf_assign_int
        module procedure cx_assign_mf
        module procedure cx_assign_dp
        module procedure cx_assign_cdp
    end interface
    public :: assignment(=)

    ! Defined I/O
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

    ! Named Constants
    type(float64x2), parameter, public :: MF_ZERO = float64x2([0.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_ONE = float64x2([1.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_TWO = float64x2([2.0_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_HALF = float64x2([0.5_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_EIGHT = float64x2([8.0_dp, 0.0_dp])
    
    ! Scaling constants
    type(float64x2), parameter, public :: MF_SAFMIN = float64x2([tiny(0.0_dp), 0.0_dp])
    type(float64x2), parameter, public :: MF_SAFMAX = float64x2([huge(0.0_dp), 0.0_dp])
    type(float64x2), parameter, public :: MF_TSML = float64x2([1.0e-100_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_TBIG = float64x2([1.0e+100_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_SSML = float64x2([1.0e-50_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_SBIG = float64x2([1.0e+50_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_RTMIN = float64x2([1.0e-150_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_RTMAX = float64x2([1.0e+150_dp, 0.0_dp])

contains

    ! Internal conversion helpers
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

    ! Constructors
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
        res = mf_from_dp(real(x, dp))
    end function

    pure function cx_from_mf(r, i) result(res)
        type(float64x2), intent(in) :: r, i
        type(complex128x2) :: res
        res%re = r
        res%im = i
    end function

    pure function cx_from_dp(r, i) result(res)
        real(dp), intent(in) :: r, i
        type(complex128x2) :: res
        res%re = from_qp(real(r, qp))
        res%im = from_qp(real(i, qp))
    end function

    pure function cx_from_mf_1(r) result(res)
        type(float64x2), intent(in) :: r
        type(complex128x2) :: res
        res%re = r
        res%im = mf_from_dp(0.0_dp)
    end function
    
    pure function cx_from_dp_1(r) result(res)
        real(dp), intent(in) :: r
        type(complex128x2) :: res
        res%re = mf_from_dp(r)
        res%im = mf_from_dp(0.0_dp)
    end function

    ! Assignment
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

    pure subroutine cx_assign_mf(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        type(float64x2), intent(in) :: rhs
        lhs%re = rhs
        lhs%im = mf_from_dp(0.0_dp)
    end subroutine

    pure subroutine mf_assign_int(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        integer, intent(in) :: rhs
        lhs = from_qp(real(rhs, qp))
    end subroutine

    pure subroutine cx_assign_dp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        real(dp), intent(in) :: rhs
        lhs%re = from_qp(real(rhs, qp))
        lhs%im = from_qp(0.0_qp)
    end subroutine

    pure subroutine cx_assign_cdp(lhs, rhs)
        type(complex128x2), intent(out) :: lhs
        complex(dp), intent(in) :: rhs
        lhs%re = from_qp(real(real(rhs, dp), qp))
        lhs%im = from_qp(real(aimag(rhs), qp))
    end subroutine

    ! Conversions
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

    ! Unary
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

    ! Math Generics (Real)
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

    pure function mf_atan2(y, x) result(res)
        type(float64x2), intent(in) :: y, x
        type(float64x2) :: res
        res = from_qp(atan2(to_qp(y), to_qp(x)))
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

    ! Math Generics (Complex)
    pure function cx_abs(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(abs(to_cqp(x)))
    end function

    pure function cx_sqrt(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(sqrt(to_cqp(x)))
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

    pure function cx_aimag(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(aimag(to_cqp(x)))
    end function

    pure function cx_conjg(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res = from_cqp(conjg(to_cqp(x)))
    end function

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
        res%im = mf_from_dp(0.0_dp)
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

    pure function mf_min(x, y) result(res)
        type(float64x2), intent(in) :: x, y
        type(float64x2) :: res
        res = from_qp(min(to_qp(x), to_qp(y)))
    end function

    pure function mf_min3(x, y, z) result(res)
        type(float64x2), intent(in) :: x, y, z
        type(float64x2) :: res
        res = from_qp(min(to_qp(x), to_qp(y), to_qp(z)))
    end function

    pure function mf_max(x, y) result(res)
        type(float64x2), intent(in) :: x, y
        type(float64x2) :: res
        res = from_qp(max(to_qp(x), to_qp(y)))
    end function

    pure function mf_max3(x, y, z) result(res)
        type(float64x2), intent(in) :: x, y, z
        type(float64x2) :: res
        res = from_qp(max(to_qp(x), to_qp(y), to_qp(z)))
    end function

    pure function mf_sign(x, y) result(res)
        type(float64x2), intent(in) :: x, y
        type(float64x2) :: res
        res = from_qp(sign(to_qp(x), to_qp(y)))
    end function

    pure function mf_mod(x, y) result(res)
        type(float64x2), intent(in) :: x, y
        type(float64x2) :: res
        res = from_qp(mod(to_qp(x), to_qp(y)))
    end function

    ! Binary Arithmetic
    pure function mf_add_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + to_qp(b))
    end function
    pure function mf_add_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + real(b, qp))
    end function
    pure function dp_add_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function
    pure function mf_add_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + real(b, qp))
    end function
    pure function sp_add_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function
    pure function mf_add_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) + b)
    end function
    pure function int_add_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) + to_qp(b))
    end function

    pure function cx_add_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a, b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + to_cqp(b))
    end function
    pure function cx_add_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) + to_qp(b))
    end function
    pure function mf_add_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) + to_cqp(b))
    end function

    pure function mf_sub_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - to_qp(b))
    end function
    pure function mf_sub_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - real(b, qp))
    end function
    pure function dp_sub_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function
    pure function mf_sub_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - real(b, qp))
    end function
    pure function sp_sub_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function
    pure function mf_sub_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) - b)
    end function
    pure function int_sub_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) - to_qp(b))
    end function

    pure function cx_sub_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a, b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - to_cqp(b))
    end function
    pure function cx_sub_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) - to_qp(b))
    end function
    pure function mf_sub_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) - to_cqp(b))
    end function

    pure function mf_mul_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * to_qp(b))
    end function
    pure function mf_mul_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * real(b, qp))
    end function
    pure function dp_mul_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function
    pure function mf_mul_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * real(b, qp))
    end function
    pure function sp_mul_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function
    pure function mf_mul_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) * b)
    end function
    pure function int_mul_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) * to_qp(b))
    end function

    pure function cx_mul_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a, b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * to_cqp(b))
    end function
    pure function cx_mul_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) * to_qp(b))
    end function
    pure function mf_mul_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) * to_cqp(b))
    end function

    pure function mf_div_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / to_qp(b))
    end function
    pure function mf_div_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / real(b, qp))
    end function
    pure function dp_div_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function
    pure function mf_div_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / real(b, qp))
    end function
    pure function sp_div_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function
    pure function mf_div_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) / b)
    end function
    pure function int_div_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res = from_qp(real(a, qp) / to_qp(b))
    end function

    pure function cx_div_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a, b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / to_cqp(b))
    end function
    pure function cx_div_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_cqp(a) / to_qp(b))
    end function
    pure function mf_div_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        type(complex128x2) :: res
        res = from_cqp(to_qp(a) / to_cqp(b))
    end function

    pure function mf_pow_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** b)
    end function

    pure function mf_pow_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** to_qp(b))
    end function

    ! Comparison
    pure function mf_eq_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = to_qp(a) == to_qp(b)
    end function
    pure function mf_eq_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function
    pure function dp_eq_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function
    pure function mf_eq_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function
    pure function sp_eq_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function
    pure function mf_eq_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) == real(b, qp)
    end function
    pure function int_eq_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == to_qp(b)
    end function
    pure function cx_eq_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a, b
        logical :: res
        res = to_cqp(a) == to_cqp(b)
    end function
    pure function cx_eq_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) == to_qp(b)
    end function
    pure function mf_eq_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_qp(a) == to_cqp(b)
    end function

    pure function mf_ne_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = to_qp(a) /= to_qp(b)
    end function
    pure function mf_ne_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function
    pure function dp_ne_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function
    pure function mf_ne_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function
    pure function sp_ne_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function
    pure function mf_ne_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) /= real(b, qp)
    end function
    pure function int_ne_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= to_qp(b)
    end function
    pure function cx_ne_cx(a, b) result(res)
        type(complex128x2), intent(in) :: a, b
        logical :: res
        res = to_cqp(a) /= to_cqp(b)
    end function
    pure function cx_ne_mf(a, b) result(res)
        type(complex128x2), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = to_cqp(a) /= to_qp(b)
    end function
    pure function mf_ne_cx(a, b) result(res)
        type(float64x2), intent(in) :: a
        type(complex128x2), intent(in) :: b
        logical :: res
        res = to_qp(a) /= to_cqp(b)
    end function

    pure function mf_lt_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = to_qp(a) < to_qp(b)
    end function
    pure function mf_lt_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function
    pure function dp_lt_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function
    pure function mf_lt_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function
    pure function sp_lt_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function
    pure function mf_lt_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) < real(b, qp)
    end function
    pure function int_lt_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < to_qp(b)
    end function

    pure function mf_gt_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = to_qp(a) > to_qp(b)
    end function
    pure function mf_gt_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function
    pure function dp_gt_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function
    pure function mf_gt_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function
    pure function sp_gt_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function
    pure function mf_gt_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) > real(b, qp)
    end function
    pure function int_gt_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > to_qp(b)
    end function

    pure function mf_le_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = to_qp(a) <= to_qp(b)
    end function
    pure function mf_le_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function
    pure function dp_le_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function
    pure function mf_le_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function
    pure function sp_le_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function
    pure function mf_le_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) <= real(b, qp)
    end function
    pure function int_le_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= to_qp(b)
    end function

    pure function mf_ge_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = to_qp(a) >= to_qp(b)
    end function
    pure function mf_ge_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function
    pure function dp_ge_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function
    pure function mf_ge_sp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(sp), intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function
    pure function sp_ge_mf(a, b) result(res)
        real(sp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function
    pure function mf_ge_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        logical :: res
        res = to_qp(a) >= real(b, qp)
    end function
    pure function int_ge_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= to_qp(b)
    end function

    ! Defined I/O Implementations
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
