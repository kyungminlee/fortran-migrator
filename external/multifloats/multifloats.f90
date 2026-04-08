module multifloats
    implicit none
    private

    integer, parameter :: qp = 16  ! Quad precision
    integer, parameter :: dp = 8   ! Double precision

    type, public :: float64x2
        real(dp) :: limbs(2)
    end type float64x2

    type, public :: complex128x2
        type(float64x2) :: re, im
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
        module procedure cx_add_cx
    end interface
    public :: operator(+)

    interface operator(-)
        module procedure mf_sub_mf
        module procedure mf_sub_dp
        module procedure dp_sub_mf
        module procedure mf_sub_int
        module procedure int_sub_mf
        module procedure mf_neg
        module procedure cx_sub_cx
        module procedure cx_neg
    end interface
    public :: operator(-)

    interface operator(*)
        module procedure mf_mul_mf
        module procedure mf_mul_dp
        module procedure dp_mul_mf
        module procedure mf_mul_int
        module procedure int_mul_mf
        module procedure cx_mul_cx
    end interface
    public :: operator(*)

    interface operator(/)
        module procedure mf_div_mf
        module procedure mf_div_dp
        module procedure dp_div_mf
        module procedure mf_div_int
        module procedure int_div_mf
        module procedure cx_div_cx
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
    
    ! Scaling constants - using largest/smallest representable DP values for mock safety
    type(float64x2), parameter, public :: MF_SAFMIN = float64x2([tiny(0.0_dp), 0.0_dp])
    type(float64x2), parameter, public :: MF_SAFMAX = float64x2([huge(0.0_dp), 0.0_dp])
    type(float64x2), parameter, public :: MF_TSML = float64x2([1.0e-100_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_TBIG = float64x2([1.0e+100_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_SSML = float64x2([1.0e-50_dp, 0.0_dp])
    type(float64x2), parameter, public :: MF_SBIG = float64x2([1.0e+50_dp, 0.0_dp])

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
    
    pure function mf_from_int(v) result(res)
        integer, intent(in) :: v
        type(float64x2) :: res
        res = from_qp(real(v, qp))
    end function

    pure function mf_from_char(v) result(res)
        character(len=*), intent(in) :: v
        type(float64x2) :: res
        real(qp) :: val_qp
        read(v, *) val_qp
        res = from_qp(val_qp)
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

    ! Assignment from DP
    pure subroutine mf_assign_dp(lhs, rhs)
        type(float64x2), intent(out) :: lhs
        real(dp), intent(in) :: rhs
        lhs = from_qp(real(rhs, qp))
    end subroutine

    ! Conversions
    pure function mf_dble(x) result(res)
        type(float64x2), intent(in) :: x
        real(dp) :: res
        res = real(to_qp(x), dp)
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

    pure function cx_real(x) result(res)
        type(complex128x2), intent(in) :: x
        type(float64x2) :: res
        res = from_qp(real(to_cqp(x), qp))
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

    pure function mf_pow_int(a, b) result(res)
        type(float64x2), intent(in) :: a
        integer, intent(in) :: b
        type(float64x2) :: res
        res = from_qp(to_qp(a) ** b)
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
            ! Use the part after DT as the format
            write(unit, fmt='(' // iotype(3:) // ')', iostat=iostat, iomsg=iomsg) val_qp
        else
            ! Default behavior for DT or other
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
        read(unit, *, iostat=iostat, iomsg=iomsg) val_qp
        dtv%limbs(1) = real(val_qp, dp)
        dtv%limbs(2) = real(val_qp - real(dtv%limbs(1), qp), dp)
    end subroutine

    subroutine read_cx_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
        class(complex128x2), intent(inout) :: dtv
        integer, intent(in) :: unit
        character(*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(*), intent(inout) :: iomsg
        complex(qp) :: val_cqp
        read(unit, *, iostat=iostat, iomsg=iomsg) val_cqp
        dtv%re%limbs(1) = real(real(val_cqp, qp), dp)
        dtv%re%limbs(2) = real(real(val_cqp, qp) - real(dtv%re%limbs(1), qp), dp)
        dtv%im%limbs(1) = real(aimag(val_cqp), dp)
        dtv%im%limbs(2) = real(aimag(val_cqp) - real(dtv%im%limbs(1), qp), dp)
    end subroutine

end module multifloats
