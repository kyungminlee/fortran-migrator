module multifloats
    implicit none
    private

    integer, parameter :: qp = 16  ! Quad precision
    integer, parameter :: dp = 8   ! Double precision

    type, public :: float64x2
        real(qp) :: val
    end type float64x2

    type, public :: complex128x2
        complex(qp) :: val
    end type complex128x2

    ! Constructors
    interface float64x2
        module procedure mf_from_dp
        module procedure mf_from_int
        module procedure mf_from_char
    end interface float64x2

    interface complex128x2
        module procedure cx_from_mf
        module procedure cx_from_dp
    end interface complex128x2

    ! Arithmetic Operators
    interface operator(+)
        module procedure mf_add_mf
        module procedure mf_add_dp
        module procedure dp_add_mf
        module procedure mf_add_int
        module procedure int_add_mf
    end interface
    public :: operator(+)

    interface operator(-)
        module procedure mf_sub_mf
        module procedure mf_sub_dp
        module procedure dp_sub_mf
        module procedure mf_sub_int
        module procedure int_sub_mf
        module procedure mf_neg
    end interface
    public :: operator(-)

    interface operator(*)
        module procedure mf_mul_mf
        module procedure mf_mul_dp
        module procedure dp_mul_mf
        module procedure mf_mul_int
        module procedure int_mul_mf
    end interface
    public :: operator(*)

    interface operator(/)
        module procedure mf_div_mf
        module procedure mf_div_dp
        module procedure dp_div_mf
        module procedure mf_div_int
        module procedure int_div_mf
    end interface
    public :: operator(/)

    interface operator(**)
        module procedure mf_pow_int
    end interface
    public :: operator(**)

    ! Comparison Operators
    interface operator(==)
        module procedure mf_eq_mf
        module procedure mf_eq_dp
        module procedure dp_eq_mf
    end interface
    public :: operator(==)

    interface operator(/=)
        module procedure mf_ne_mf
        module procedure mf_ne_dp
        module procedure dp_ne_mf
    end interface
    public :: operator(/=)

    interface operator(<)
        module procedure mf_lt_mf
        module procedure mf_lt_dp
        module procedure dp_lt_mf
    end interface
    public :: operator(<)

    interface operator(>)
        module procedure mf_gt_mf
        module procedure mf_gt_dp
        module procedure dp_gt_mf
    end interface
    public :: operator(>)

    interface operator(<=)
        module procedure mf_le_mf
        module procedure mf_le_dp
        module procedure dp_le_mf
    end interface
    public :: operator(<=)

    interface operator(>=)
        module procedure mf_ge_mf
        module procedure mf_ge_dp
        module procedure dp_ge_mf
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
    end interface
    public :: sin

    interface cos
        module procedure mf_cos
    end interface
    public :: cos

    interface tan
        module procedure mf_tan
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

    interface aimag
        module procedure cx_aimag
    end interface
    public :: aimag

    interface conjg
        module procedure cx_conjg
    end interface
    public :: conjg

    interface real
        module procedure cx_real
    end interface
    public :: real

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

    public :: mf_to_double

    interface assignment(=)
        module procedure mf_assign_dp
    end interface
    public :: assignment(=)

    ! Named Constants
    type(float64x2), parameter, public :: MF_ZERO = float64x2(0.0_qp)
    type(float64x2), parameter, public :: MF_ONE = float64x2(1.0_qp)
    type(float64x2), parameter, public :: MF_TWO = float64x2(2.0_qp)
    type(float64x2), parameter, public :: MF_HALF = float64x2(0.5_qp)
    type(float64x2), parameter, public :: MF_EIGHT = float64x2(8.0_qp)
    
    ! Dummy scaling constants for mock (quad precision approximations)
    type(float64x2), parameter, public :: MF_SAFMIN = float64x2(tiny(0.0_qp))
    type(float64x2), parameter, public :: MF_SAFMAX = float64x2(huge(0.0_qp))
    type(float64x2), parameter, public :: MF_TSML = float64x2(1.0e-100_qp)
    type(float64x2), parameter, public :: MF_TBIG = float64x2(1.0e+100_qp)
    type(float64x2), parameter, public :: MF_SSML = float64x2(1.0e-50_qp)
    type(float64x2), parameter, public :: MF_SBIG = float64x2(1.0e+50_qp)

contains

    ! Constructors
    pure function mf_from_dp(v) result(res)
        real(dp), intent(in) :: v
        type(float64x2) :: res
        res%val = real(v, qp)
    end function
    
    pure function mf_from_int(v) result(res)
        integer, intent(in) :: v
        type(float64x2) :: res
        res%val = real(v, qp)
    end function

    pure function mf_from_char(v) result(res)
        character(len=*), intent(in) :: v
        type(float64x2) :: res
        read(v, *) res%val
    end function

    pure function cx_from_mf(r, i) result(res)
        type(float64x2), intent(in) :: r, i
        type(complex128x2) :: res
        res%val = cmplx(r%val, i%val, kind=qp)
    end function

    pure function cx_from_dp(r, i) result(res)
        real(dp), intent(in) :: r, i
        type(complex128x2) :: res
        res%val = cmplx(real(r, qp), real(i, qp), kind=qp)
    end function

    ! Assignment
    pure subroutine mf_assign_dp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        real(dp), intent(in) :: rhs
        lhs%val = real(rhs, qp)
    end subroutine

    ! Conversions
    pure function mf_dble(x) result(res)
        type(float64x2), intent(in) :: x
        real(dp) :: res
        res = real(x%val, dp)
    end function

    pure function mf_to_double(x) result(res)
        type(float64x2), intent(in) :: x
        real(dp) :: res
        res = real(x%val, dp)
    end function

    ! Unary
    pure function mf_neg(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = -x%val
    end function

    ! Math Generics (Real)
    pure function mf_abs(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = abs(x%val)
    end function

    pure function mf_sqrt(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = sqrt(x%val)
    end function

    pure function mf_sin(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = sin(x%val)
    end function

    pure function mf_cos(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = cos(x%val)
    end function

    pure function mf_tan(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = tan(x%val)
    end function

    pure function mf_exp(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = exp(x%val)
    end function

    pure function mf_log(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = log(x%val)
    end function

    pure function mf_log10(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res
        res%val = log10(x%val)
    end function

    ! Math Generics (Complex)
    pure function cx_abs(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res%val = abs(x%val)
    end function

    pure function cx_sqrt(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res%val = sqrt(x%val)
    end function

    pure function cx_exp(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res%val = exp(x%val)
    end function

    pure function cx_log(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res%val = log(x%val)
    end function

    pure function cx_aimag(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res%val = aimag(x%val)
    end function

    pure function cx_conjg(x) result(res)
        type(complex128x2), intent(in) :: x
        type(complex128x2) :: res
        res%val = conjg(x%val)
    end function

    pure function cx_real(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res%val = real(x%val, qp)
    end function

    pure function mf_min(x, y) result(res)
        type(float64x2), intent(in) :: x, y
        type(float64x2) :: res
        res%val = min(x%val, y%val)
    end function

    pure function mf_min3(x, y, z) result(res)
        type(float64x2), intent(in) :: x, y, z
        type(float64x2) :: res
        res%val = min(x%val, y%val, z%val)
    end function

    pure function mf_max(x, y) result(res)
        type(float64x2), intent(in) :: x, y
        type(float64x2) :: res
        res%val = max(x%val, y%val)
    end function

    pure function mf_max3(x, y, z) result(res)
        type(float64x2), intent(in) :: x, y, z
        type(float64x2) :: res
        res%val = max(x%val, y%val, z%val)
    end function

    pure function mf_sign(x, y) result(res)
        type(float64x2), intent(in) :: x, y
        type(float64x2) :: res
        res%val = sign(x%val, y%val)
    end function

    pure function mf_mod(x, y) result(res)
        type(float64x2), intent(in) :: x, y
        type(float64x2) :: res
        res%val = mod(x%val, y%val)
    end function

    ! Binary Arithmetic (Macro-style functions)
    ! Add
    pure function mf_add_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res%val = a%val + b%val
    end function
    pure function mf_add_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res%val = a%val + real(b, qp)
    end function
    pure function dp_add_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res%val = real(a, qp) + b%val
    end function
    pure function mf_add_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res%val = a%val + b
    end function
    pure function int_add_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res%val = a + b%val
    end function

    ! Sub
    pure function mf_sub_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res%val = a%val - b%val
    end function
    pure function mf_sub_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res%val = a%val - real(b, qp)
    end function
    pure function dp_sub_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res%val = real(a, qp) - b%val
    end function
    pure function mf_sub_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res%val = a%val - b
    end function
    pure function int_sub_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res%val = a - b%val
    end function

    ! Mul
    pure function mf_mul_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res%val = a%val * b%val
    end function
    pure function mf_mul_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res%val = a%val * real(b, qp)
    end function
    pure function dp_mul_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res%val = real(a, qp) * b%val
    end function
    pure function mf_mul_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res%val = a%val * b
    end function
    pure function int_mul_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res%val = a * b%val
    end function

    ! Div
    pure function mf_div_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        type(float64x2) :: res
        res%val = a%val / b%val
    end function
    pure function mf_div_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        type(float64x2) :: res
        res%val = a%val / real(b, qp)
    end function
    pure function dp_div_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res%val = real(a, qp) / b%val
    end function
    pure function mf_div_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res%val = a%val / b
    end function
    pure function int_div_mf(a, b) result(res)
        integer, intent(in) :: a
        type(float64x2), intent(in) :: b
        type(float64x2) :: res
        res%val = a / b%val
    end function

    ! Pow
    pure function mf_pow_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res%val = a%val ** b
    end function

    ! Comparison
    pure function mf_eq_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = a%val == b%val
    end function
    pure function mf_eq_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = a%val == real(b, qp)
    end function
    pure function dp_eq_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) == b%val
    end function

    pure function mf_ne_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = a%val /= b%val
    end function
    pure function mf_ne_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = a%val /= real(b, qp)
    end function
    pure function dp_ne_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) /= b%val
    end function

    pure function mf_lt_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = a%val < b%val
    end function
    pure function mf_lt_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = a%val < real(b, qp)
    end function
    pure function dp_lt_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) < b%val
    end function

    pure function mf_gt_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = a%val > b%val
    end function
    pure function mf_gt_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = a%val > real(b, qp)
    end function
    pure function dp_gt_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) > b%val
    end function

    pure function mf_le_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = a%val <= b%val
    end function
    pure function mf_le_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = a%val <= real(b, qp)
    end function
    pure function dp_le_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) <= b%val
    end function

    pure function mf_ge_mf(a, b) result(res)
        type(float64x2), intent(in) :: a, b
        logical :: res
        res = a%val >= b%val
    end function
    pure function mf_ge_dp(a, b) result(res)
        type(float64x2), intent(in) :: a
        real(dp), intent(in) :: b
        logical :: res
        res = a%val >= real(b, qp)
    end function
    pure function dp_ge_mf(a, b) result(res)
        real(dp), intent(in) :: a
        type(float64x2), intent(in) :: b
        logical :: res
        res = real(a, qp) >= b%val
    end function

end module multifloats
