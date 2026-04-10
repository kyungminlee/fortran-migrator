program multifloats_unittest
    use multifloats
    implicit none

    integer :: num_tests = 0
    integer :: num_passed = 0

    print *, "Starting Multifloats Unit Tests..."
    print *, "Internal structure: real(8) limbs(2), Op precision: real(16)"
    print *, "--------------------------------------------------------"

    call test_constructors()
    call test_assignment()
    call test_arithmetic_mf_mf()
    call test_arithmetic_mf_dp()
    call test_arithmetic_mf_sp()
    call test_arithmetic_mf_int()
    call test_arithmetic_complex()
    call test_arithmetic_mixed_complex()
    call test_power()
    call test_comparisons_mf()
    call test_comparisons_int()
    call test_comparisons_complex()
    call test_math_generics()
    call test_math_trig_inverse()
    call test_binary_real_funcs()
    call test_binary_real_mixed()
    call test_complex_funcs()
    call test_conversions()
    call test_io()

    print *, "--------------------------------------------------------"
    print *, "Tests Run:   ", num_tests
    print *, "Tests Passed:", num_passed
    if (num_tests == num_passed) then
        print *, "RESULT: ALL TESTS PASSED"
    else
        print *, "RESULT: SOME TESTS FAILED"
        stop 1
    end if

contains

    subroutine check(condition, msg)
        logical, intent(in) :: condition
        character(*), intent(in) :: msg
        num_tests = num_tests + 1
        if (condition) then
            num_passed = num_passed + 1
        else
            print *, "[FAIL] ", msg
        end if
    end subroutine

    ! ================================================================
    ! Constructors
    ! ================================================================

    subroutine test_constructors()
        type(float64x2) :: a
        real(8) :: dval = 1.2345678901234567_8
        real(4) :: sval = 1.5
        complex(8) :: cdpval = (3.0_8, 4.0_8)
        complex(4) :: cspval = (1.5, 2.5)

        a = float64x2(dval)
        call check(abs(mf_to_double(a) - dval) < 1e-15, "Ctor: from DP")

        a = float64x2(sval)
        call check(abs(mf_to_double(a) - 1.5_8) < 1e-6, "Ctor: from SP")

        a = float64x2(42)
        call check(abs(mf_to_double(a) - 42.0_8) < 1e-15, "Ctor: from INT")

        a = float64x2('3.14159265358979323846')
        call check(abs(mf_to_double(a) - 3.14159265358979_8) < 1e-14, "Ctor: from CHAR")

        ! From complex(dp) — takes real part
        a = float64x2(cdpval)
        call check(abs(mf_to_double(a) - 3.0_8) < 1e-15, "Ctor: from complex(dp)")

        ! From complex(sp) — takes real part
        a = float64x2(cspval)
        call check(abs(mf_to_double(a) - 1.5_8) < 1e-6, "Ctor: from complex(sp)")

        ! From complex128x2 — takes real part
        block
            type(complex128x2) :: z
            z = complex128x2(float64x2(7.0_8), float64x2(9.0_8))
            a = float64x2(z)
            call check(mf_to_double(a) == 7.0_8, "Ctor: from complex128x2")
        end block
    end subroutine

    ! ================================================================
    ! Assignment
    ! ================================================================

    subroutine test_assignment()
        type(float64x2) :: m
        type(complex128x2) :: z

        ! mf = dp
        m = 3.5_8
        call check(mf_to_double(m) == 3.5_8, "Assign: mf = dp")

        ! mf = sp
        m = 2.5
        call check(abs(mf_to_double(m) - 2.5_8) < 1e-6, "Assign: mf = sp")

        ! mf = int
        m = 7
        call check(mf_to_double(m) == 7.0_8, "Assign: mf = int")

        ! mf = complex128x2 (takes real part)
        z = complex128x2(float64x2(11.0_8), float64x2(22.0_8))
        m = z
        call check(mf_to_double(m) == 11.0_8, "Assign: mf = cx (real part)")

        ! mf = complex(dp) (takes real part)
        m = (5.0_8, 6.0_8)
        call check(mf_to_double(m) == 5.0_8, "Assign: mf = cdp (real part)")

        ! mf = complex(sp) (takes real part)
        m = cmplx(1.5, 2.5)
        call check(abs(mf_to_double(m) - 1.5_8) < 1e-6, "Assign: mf = csp (real part)")

        ! cx = mf
        z = float64x2(4.0_8)
        call check(mf_to_double(real(z)) == 4.0_8, "Assign: cx = mf")
        call check(mf_to_double(aimag(z)) == 0.0_8, "Assign: cx = mf (imag=0)")

        ! cx = dp
        z = 9.0_8
        call check(mf_to_double(real(z)) == 9.0_8, "Assign: cx = dp")

        ! cx = sp
        z = 3.0
        call check(abs(mf_to_double(real(z)) - 3.0_8) < 1e-6, "Assign: cx = sp")

        ! cx = int
        z = 5
        call check(mf_to_double(real(z)) == 5.0_8, "Assign: cx = int")

        ! cx = complex(dp)
        z = (2.0_8, 3.0_8)
        call check(mf_to_double(real(z)) == 2.0_8, "Assign: cx = cdp (re)")
        call check(mf_to_double(aimag(z)) == 3.0_8, "Assign: cx = cdp (im)")

        ! cx = complex(sp)
        z = cmplx(1.0, 2.0)
        call check(abs(mf_to_double(real(z)) - 1.0_8) < 1e-6, "Assign: cx = csp (re)")
        call check(abs(mf_to_double(aimag(z)) - 2.0_8) < 1e-6, "Assign: cx = csp (im)")
    end subroutine

    ! ================================================================
    ! Arithmetic: mf OP mf
    ! ================================================================

    subroutine test_arithmetic_mf_mf()
        type(float64x2) :: a, b, r

        a = float64x2(10.0_8)
        b = float64x2(5.0_8)

        r = a + b;  call check(mf_to_double(r) == 15.0_8, "Arith: mf + mf")
        r = a - b;  call check(mf_to_double(r) == 5.0_8,  "Arith: mf - mf")
        r = a * b;  call check(mf_to_double(r) == 50.0_8, "Arith: mf * mf")
        r = a / b;  call check(mf_to_double(r) == 2.0_8,  "Arith: mf / mf")
        r = -a;     call check(mf_to_double(r) == -10.0_8, "Arith: -mf")
    end subroutine

    ! ================================================================
    ! Arithmetic: mf OP dp, dp OP mf
    ! ================================================================

    subroutine test_arithmetic_mf_dp()
        type(float64x2) :: a, r
        real(8) :: d = 2.0_8

        a = float64x2(10.0_8)

        r = a + d;  call check(mf_to_double(r) == 12.0_8, "Arith: mf + dp")
        r = a - d;  call check(mf_to_double(r) == 8.0_8,  "Arith: mf - dp")
        r = a * d;  call check(mf_to_double(r) == 20.0_8, "Arith: mf * dp")
        r = a / d;  call check(mf_to_double(r) == 5.0_8,  "Arith: mf / dp")

        r = d + a;  call check(mf_to_double(r) == 12.0_8, "Arith: dp + mf")
        r = d - a;  call check(mf_to_double(r) == -8.0_8, "Arith: dp - mf")
        r = d * a;  call check(mf_to_double(r) == 20.0_8, "Arith: dp * mf")
        r = d / a;  call check(mf_to_double(r) == 0.2_8,  "Arith: dp / mf")
    end subroutine

    ! ================================================================
    ! Arithmetic: mf OP sp, sp OP mf
    ! ================================================================

    subroutine test_arithmetic_mf_sp()
        type(float64x2) :: a, r
        real(4) :: s = 3.0

        a = float64x2(10.0_8)

        r = a + s;  call check(mf_to_double(r) == 13.0_8, "Arith: mf + sp")
        r = a - s;  call check(mf_to_double(r) == 7.0_8,  "Arith: mf - sp")
        r = a * s;  call check(mf_to_double(r) == 30.0_8, "Arith: mf * sp")
        r = a / s;  call check(abs(mf_to_double(r) - 3.333333333333_8) < 1e-10, "Arith: mf / sp")

        r = s + a;  call check(mf_to_double(r) == 13.0_8, "Arith: sp + mf")
        r = s - a;  call check(mf_to_double(r) == -7.0_8, "Arith: sp - mf")
        r = s * a;  call check(mf_to_double(r) == 30.0_8, "Arith: sp * mf")
    end subroutine

    ! ================================================================
    ! Arithmetic: mf OP int, int OP mf
    ! ================================================================

    subroutine test_arithmetic_mf_int()
        type(float64x2) :: a, r
        integer :: i = 3

        a = float64x2(10.0_8)

        r = a + i;  call check(mf_to_double(r) == 13.0_8, "Arith: mf + int")
        r = a - i;  call check(mf_to_double(r) == 7.0_8,  "Arith: mf - int")
        r = a * i;  call check(mf_to_double(r) == 30.0_8, "Arith: mf * int")
        r = a / i;  call check(abs(mf_to_double(r) - 3.33333333333333_8) < 1e-14, "Arith: mf / int")

        r = i + a;  call check(mf_to_double(r) == 13.0_8, "Arith: int + mf")
        r = i - a;  call check(mf_to_double(r) == -7.0_8, "Arith: int - mf")
        r = i * a;  call check(mf_to_double(r) == 30.0_8, "Arith: int * mf")
        r = i / a;  call check(mf_to_double(r) == 0.3_8,  "Arith: int / mf (approx)")
    end subroutine

    ! ================================================================
    ! Arithmetic: cx OP cx, cx OP mf, mf OP cx
    ! ================================================================

    subroutine test_arithmetic_complex()
        type(float64x2) :: r_part
        type(complex128x2) :: z1, z2, zr

        z1 = complex128x2(float64x2(3.0_8), float64x2(4.0_8))
        z2 = complex128x2(float64x2(1.0_8), float64x2(2.0_8))

        ! cx + cx
        zr = z1 + z2
        call check(mf_to_double(real(zr)) == 4.0_8, "Arith: cx + cx (re)")
        call check(mf_to_double(aimag(zr)) == 6.0_8, "Arith: cx + cx (im)")

        ! cx - cx
        zr = z1 - z2
        call check(mf_to_double(real(zr)) == 2.0_8, "Arith: cx - cx (re)")
        call check(mf_to_double(aimag(zr)) == 2.0_8, "Arith: cx - cx (im)")

        ! cx * cx: (3+4i)*(1+2i) = (3-8) + (6+4)i = -5+10i
        zr = z1 * z2
        call check(mf_to_double(real(zr)) == -5.0_8, "Arith: cx * cx (re)")
        call check(mf_to_double(aimag(zr)) == 10.0_8, "Arith: cx * cx (im)")

        ! cx + mf
        r_part = float64x2(10.0_8)
        zr = z1 + r_part
        call check(mf_to_double(real(zr)) == 13.0_8, "Arith: cx + mf (re)")
        call check(mf_to_double(aimag(zr)) == 4.0_8,  "Arith: cx + mf (im)")

        ! mf + cx
        zr = r_part + z1
        call check(mf_to_double(real(zr)) == 13.0_8, "Arith: mf + cx (re)")

        ! -cx
        zr = -z1
        call check(mf_to_double(real(zr)) == -3.0_8, "Arith: -cx (re)")
        call check(mf_to_double(aimag(zr)) == -4.0_8, "Arith: -cx (im)")
    end subroutine

    ! ================================================================
    ! Arithmetic: mixed cx with dp, sp, int, cdp, csp
    ! ================================================================

    subroutine test_arithmetic_mixed_complex()
        type(complex128x2) :: z, zr
        real(8) :: d = 2.0_8
        real(4) :: s = 3.0
        integer :: i = 5
        complex(8) :: cdpv = (1.0_8, 2.0_8)
        complex(4) :: cspv = (1.0, 1.0)

        z = complex128x2(float64x2(4.0_8), float64x2(3.0_8))

        ! cx + dp / dp + cx
        zr = z + d
        call check(mf_to_double(real(zr)) == 6.0_8, "Mixed: cx + dp")
        zr = d + z
        call check(mf_to_double(real(zr)) == 6.0_8, "Mixed: dp + cx")

        ! cx * sp
        zr = z * s
        call check(mf_to_double(real(zr)) == 12.0_8, "Mixed: cx * sp (re)")
        call check(mf_to_double(aimag(zr)) == 9.0_8,  "Mixed: cx * sp (im)")

        ! cx + int
        zr = z + i
        call check(mf_to_double(real(zr)) == 9.0_8, "Mixed: cx + int")

        ! int * cx
        zr = i * z
        call check(mf_to_double(real(zr)) == 20.0_8, "Mixed: int * cx (re)")
        call check(mf_to_double(aimag(zr)) == 15.0_8, "Mixed: int * cx (im)")

        ! cx + cdp: (4+3i) + (1+2i) = (5+5i)
        zr = z + cdpv
        call check(mf_to_double(real(zr)) == 5.0_8, "Mixed: cx + cdp (re)")
        call check(mf_to_double(aimag(zr)) == 5.0_8, "Mixed: cx + cdp (im)")

        ! cdp * cx
        zr = cdpv * z
        call check(abs(mf_to_double(real(zr)) - (-2.0_8)) < 1e-14, "Mixed: cdp * cx (re)")

        ! cx + csp
        zr = z + cspv
        call check(abs(mf_to_double(real(zr)) - 5.0_8) < 1e-6, "Mixed: cx + csp (re)")

        ! mf * cx (result is complex)
        block
            type(float64x2) :: m
            m = float64x2(2.0_8)
            zr = m * z
            call check(mf_to_double(real(zr)) == 8.0_8, "Mixed: mf * cx (re)")
            call check(mf_to_double(aimag(zr)) == 6.0_8, "Mixed: mf * cx (im)")
        end block
    end subroutine

    ! ================================================================
    ! Power
    ! ================================================================

    subroutine test_power()
        type(float64x2) :: a, r

        a = float64x2(5.0_8)

        ! mf ** int (integer exponentiation)
        r = a ** 2
        call check(mf_to_double(r) == 25.0_8, "Power: mf ** int")

        r = a ** 0
        call check(mf_to_double(r) == 1.0_8, "Power: mf ** 0")

        ! mf ** mf
        r = float64x2(4.0_8) ** float64x2(0.5_8)
        call check(abs(mf_to_double(r) - 2.0_8) < 1e-14, "Power: mf ** mf (sqrt)")

        ! mf ** dp
        r = a ** 2.0_8
        call check(mf_to_double(r) == 25.0_8, "Power: mf ** dp")

        ! dp ** mf
        r = 2.0_8 ** float64x2(3.0_8)
        call check(abs(mf_to_double(r) - 8.0_8) < 1e-14, "Power: dp ** mf")

        ! int ** mf
        r = 2 ** float64x2(10.0_8)
        call check(abs(mf_to_double(r) - 1024.0_8) < 1e-10, "Power: int ** mf")

        ! cx ** int
        block
            type(complex128x2) :: z, zr
            z = complex128x2(float64x2(0.0_8), float64x2(1.0_8))  ! i
            zr = z ** 2  ! i^2 = -1
            call check(abs(mf_to_double(real(zr)) - (-1.0_8)) < 1e-14, "Power: cx ** int (re)")
            call check(abs(mf_to_double(aimag(zr))) < 1e-14, "Power: cx ** int (im)")
        end block
    end subroutine

    ! ================================================================
    ! Comparisons: mf vs mf, dp, sp
    ! ================================================================

    subroutine test_comparisons_mf()
        type(float64x2) :: a, b
        real(8) :: d = 5.0_8
        real(4) :: s = 5.0

        a = float64x2(10.0_8)
        b = float64x2(5.0_8)

        call check(a > b,  "Cmp: mf > mf")
        call check(b < a,  "Cmp: mf < mf")
        call check(a == a, "Cmp: mf == mf")
        call check(a /= b, "Cmp: mf /= mf")
        call check(a >= b, "Cmp: mf >= mf")
        call check(b <= a, "Cmp: mf <= mf")

        call check(a > d,  "Cmp: mf > dp")
        call check(b == d, "Cmp: mf == dp")
        call check(d < a,  "Cmp: dp < mf")
        call check(d == b, "Cmp: dp == mf")

        call check(a > s,  "Cmp: mf > sp")
        call check(s < a,  "Cmp: sp < mf")
        call check(s == b, "Cmp: sp == mf")
    end subroutine

    ! ================================================================
    ! Comparisons: mf vs int
    ! ================================================================

    subroutine test_comparisons_int()
        type(float64x2) :: a

        a = float64x2(10.0_8)

        call check(a == 10,  "Cmp: mf == int")
        call check(a /= 5,   "Cmp: mf /= int")
        call check(a > 5,    "Cmp: mf > int")
        call check(a < 20,   "Cmp: mf < int")
        call check(a >= 10,  "Cmp: mf >= int")
        call check(a <= 10,  "Cmp: mf <= int")

        call check(10 == a,  "Cmp: int == mf")
        call check(5 /= a,   "Cmp: int /= mf")
        call check(5 < a,    "Cmp: int < mf")
        call check(20 > a,   "Cmp: int > mf")
    end subroutine

    ! ================================================================
    ! Comparisons: complex equality
    ! ================================================================

    subroutine test_comparisons_complex()
        type(complex128x2) :: z1, z2
        type(float64x2) :: m

        z1 = complex128x2(float64x2(1.0_8), float64x2(2.0_8))
        z2 = complex128x2(float64x2(1.0_8), float64x2(2.0_8))
        m = float64x2(1.0_8)

        call check(z1 == z2, "Cmp: cx == cx")
        call check(.not. (z1 /= z2), "Cmp: cx /= cx (equal)")

        ! cx == mf: true only if imaginary is 0
        z2 = complex128x2(float64x2(1.0_8), float64x2(0.0_8))
        call check(z2 == m, "Cmp: cx == mf (im=0)")
        call check(z1 /= m, "Cmp: cx /= mf (im/=0)")

        ! cx == cdp
        call check(z2 == (1.0_8, 0.0_8), "Cmp: cx == cdp")
        call check(z1 /= (1.0_8, 0.0_8), "Cmp: cx /= cdp")
    end subroutine

    ! ================================================================
    ! Math generics: basic
    ! ================================================================

    subroutine test_math_generics()
        type(float64x2) :: a

        a = float64x2(-1.0_8)
        call check(mf_to_double(abs(a)) == 1.0_8, "Math: ABS")

        a = float64x2(4.0_8)
        call check(mf_to_double(sqrt(a)) == 2.0_8, "Math: SQRT")

        a = float64x2(0.0_8)
        call check(mf_to_double(sin(a)) == 0.0_8, "Math: SIN(0)")
        call check(mf_to_double(cos(a)) == 1.0_8, "Math: COS(0)")
        call check(mf_to_double(tan(a)) == 0.0_8, "Math: TAN(0)")

        a = float64x2(1.0_8)
        call check(abs(mf_to_double(exp(a)) - 2.71828182845904_8) < 1e-14, "Math: EXP")
        call check(mf_to_double(log(exp(a))) == 1.0_8, "Math: LOG(EXP(1))")

        a = float64x2(100.0_8)
        call check(mf_to_double(log10(a)) == 2.0_8, "Math: LOG10(100)")
    end subroutine

    ! ================================================================
    ! Math: inverse trig (atan, asin, acos, atan2)
    ! ================================================================

    subroutine test_math_trig_inverse()
        type(float64x2) :: a, r
        real(8) :: pi_8 = 3.14159265358979323846_8

        a = float64x2(0.0_8)
        call check(abs(mf_to_double(atan(a))) < 1e-15, "Math: ATAN(0)")
        call check(abs(mf_to_double(asin(a))) < 1e-15, "Math: ASIN(0)")
        call check(abs(mf_to_double(acos(a)) - pi_8/2.0_8) < 1e-14, "Math: ACOS(0)")

        a = float64x2(1.0_8)
        call check(abs(mf_to_double(atan(a)) - pi_8/4.0_8) < 1e-14, "Math: ATAN(1)")
        call check(abs(mf_to_double(asin(a)) - pi_8/2.0_8) < 1e-14, "Math: ASIN(1)")
        call check(abs(mf_to_double(acos(a))) < 1e-14, "Math: ACOS(1)")

        ! atan2(1, 1) = pi/4
        r = atan2(float64x2(1.0_8), float64x2(1.0_8))
        call check(abs(mf_to_double(r) - pi_8/4.0_8) < 1e-14, "Math: ATAN2(1,1)")

        ! atan2(0, -1) = pi
        r = atan2(float64x2(0.0_8), float64x2(-1.0_8))
        call check(abs(mf_to_double(r) - pi_8) < 1e-14, "Math: ATAN2(0,-1)")

        ! Complex atan, asin, acos (smoke test — just check they don't crash)
        block
            type(complex128x2) :: z, zr
            z = complex128x2(float64x2(0.5_8), float64x2(0.5_8))
            zr = atan(z);  call check(.true., "Math: cx ATAN (no crash)")
            zr = asin(z);  call check(.true., "Math: cx ASIN (no crash)")
            zr = acos(z);  call check(.true., "Math: cx ACOS (no crash)")
        end block
    end subroutine

    ! ================================================================
    ! Binary real functions: sign, mod, min, max
    ! ================================================================

    subroutine test_binary_real_funcs()
        type(float64x2) :: a, b, c

        a = float64x2(1.0_8)
        b = float64x2(2.0_8)
        c = float64x2(3.0_8)

        call check(mf_to_double(min(a, b)) == 1.0_8, "Math: MIN(2)")
        call check(mf_to_double(min(a, b, c)) == 1.0_8, "Math: MIN(3)")
        call check(mf_to_double(max(a, b)) == 2.0_8, "Math: MAX(2)")
        call check(mf_to_double(max(a, b, c)) == 3.0_8, "Math: MAX(3)")

        call check(mf_to_double(sign(float64x2(5.0_8), float64x2(-1.0_8))) == -5.0_8, "Math: SIGN")
        call check(mf_to_double(mod(float64x2(7.0_8), float64x2(3.0_8))) == 1.0_8, "Math: MOD")
    end subroutine

    ! ================================================================
    ! Binary real functions: mixed types
    ! ================================================================

    subroutine test_binary_real_mixed()
        type(float64x2) :: a, r
        real(8) :: d = 2.0_8
        integer :: i = 3

        a = float64x2(10.0_8)

        ! min/max: mf vs dp
        r = min(a, d)
        call check(mf_to_double(r) == 2.0_8, "Mixed: min(mf, dp)")
        r = max(d, a)
        call check(mf_to_double(r) == 10.0_8, "Mixed: max(dp, mf)")

        ! min/max: mf vs int
        r = min(a, i)
        call check(mf_to_double(r) == 3.0_8, "Mixed: min(mf, int)")
        r = max(i, a)
        call check(mf_to_double(r) == 10.0_8, "Mixed: max(int, mf)")

        ! sign: mf vs dp
        r = sign(a, -1.0_8)
        call check(mf_to_double(r) == -10.0_8, "Mixed: sign(mf, dp)")

        ! mod: mf vs dp
        r = mod(a, 3.0_8)
        call check(mf_to_double(r) == 1.0_8, "Mixed: mod(mf, dp)")

        ! atan2: mixed
        r = atan2(float64x2(1.0_8), 1.0_8)
        call check(abs(mf_to_double(r) - 0.78539816339744_8) < 1e-14, "Mixed: atan2(mf, dp)")
    end subroutine

    ! ================================================================
    ! Complex-specific functions: cmplx, real, mf_real, conjg, aimag
    ! ================================================================

    subroutine test_complex_funcs()
        type(float64x2) :: re, im, r
        type(complex128x2) :: z, zr

        re = float64x2(3.0_8)
        im = float64x2(4.0_8)

        ! cmplx(mf, mf)
        z = cmplx(re, im)
        call check(mf_to_double(real(z)) == 3.0_8, "Func: cmplx(mf,mf) re")
        call check(mf_to_double(aimag(z)) == 4.0_8, "Func: cmplx(mf,mf) im")

        ! cmplx(mf) — imaginary = 0
        z = cmplx(re)
        call check(mf_to_double(real(z)) == 3.0_8, "Func: cmplx(mf) re")
        call check(mf_to_double(aimag(z)) == 0.0_8, "Func: cmplx(mf) im=0")

        ! real(cx) and aimag(cx)
        z = complex128x2(float64x2(7.0_8), float64x2(11.0_8))
        call check(mf_to_double(real(z)) == 7.0_8, "Func: real(cx)")
        call check(mf_to_double(aimag(z)) == 11.0_8, "Func: aimag(cx)")

        ! real(mf) — identity
        r = real(re)
        call check(mf_to_double(r) == 3.0_8, "Func: real(mf) identity")

        ! mf_real(cx) — direct extraction
        r = mf_real(z)
        call check(mf_to_double(r) == 7.0_8, "Func: mf_real(cx)")

        ! mf_real(mf) — identity
        r = mf_real(re)
        call check(mf_to_double(r) == 3.0_8, "Func: mf_real(mf)")

        ! conjg
        z = complex128x2(float64x2(1.0_8), float64x2(2.0_8))
        zr = conjg(z)
        call check(mf_to_double(aimag(zr)) == -2.0_8, "Func: conjg")

        ! abs(cx) = magnitude
        z = complex128x2(float64x2(3.0_8), float64x2(4.0_8))
        call check(mf_to_double(abs(z)) == 5.0_8, "Func: abs(cx)")

        ! Complex sqrt: sqrt(i) = (1+i)/sqrt(2)
        z = complex128x2(float64x2(0.0_8), float64x2(1.0_8))
        zr = sqrt(z)
        call check(abs(mf_to_double(real(zr)) - 0.70710678118654_8) < 1e-14, "Func: cx sqrt (re)")
        call check(abs(mf_to_double(aimag(zr)) - 0.70710678118654_8) < 1e-14, "Func: cx sqrt (im)")
    end subroutine

    ! ================================================================
    ! Conversions: dble, int, nint, mf_to_double
    ! ================================================================

    subroutine test_conversions()
        type(float64x2) :: a
        real(8) :: d
        integer :: i

        a = float64x2(3.7_8)

        d = dble(a)
        call check(abs(d - 3.7_8) < 1e-15, "Conv: dble(mf)")

        d = mf_to_double(a)
        call check(abs(d - 3.7_8) < 1e-15, "Conv: mf_to_double(mf)")

        i = int(a)
        call check(i == 3, "Conv: int(mf) truncation")

        a = float64x2(3.5_8)
        i = nint(a)
        call check(i == 4, "Conv: nint(3.5)")

        a = float64x2(-2.7_8)
        i = nint(a)
        call check(i == -3, "Conv: nint(-2.7)")
    end subroutine

    ! ================================================================
    ! I/O
    ! ================================================================

    subroutine test_io()
        type(float64x2) :: a, b
        character(len=100) :: buffer

        a = float64x2(1.23456789_8)

        ! Test write with DT and specific sub-format
        write(buffer, '(DT"F12.8")') a
        call check(adjustl(trim(buffer)) == "1.23456789", "IO: WRITE (DT)")

        ! Test read
        read(buffer, '(DT)') b
        call check(abs(mf_to_double(a) - mf_to_double(b)) < 1e-8, "IO: READ (DT)")

        ! Test list-directed (smoke test)
        write(buffer, *) a
        call check(len_trim(buffer) > 0, "IO: list-directed WRITE")
    end subroutine

end program multifloats_unittest
