program multifloats_unittest
    use multifloats
    implicit none

    integer :: num_tests = 0
    integer :: num_passed = 0

    print *, "Starting Multifloats Unit Tests..."
    print *, "Internal structure: real(8) limbs(2), Op precision: real(16)"
    print *, "--------------------------------------------------------"

    call test_constructors()
    call test_arithmetic()
    call test_comparisons()
    call test_math_generics()
    call test_complex()
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
            ! print *, "[PASS] ", msg
        else
            print *, "[FAIL] ", msg
        end if
    end subroutine

    subroutine test_constructors()
        type(float64x2) :: a, b, c
        real(8) :: dval = 1.2345678901234567_8
        
        a = float64x2(dval)
        call check(abs(mf_to_double(a) - dval) < 1e-15, "Constructor from DP")

        b = float64x2(42)
        call check(abs(mf_to_double(b) - 42.0_8) < 1e-15, "Constructor from INT")

        c = float64x2('3.14159265358979323846')
        call check(abs(mf_to_double(c) - 3.14159265358979_8) < 1e-14, "Constructor from CHAR")
    end subroutine

    subroutine test_arithmetic()
        type(float64x2) :: a, b, r
        real(8) :: d = 2.0_8
        integer :: i = 3

        a = float64x2(10.0_8)
        b = float64x2(5.0_8)

        ! Real-Real
        r = a + b
        call check(mf_to_double(r) == 15.0_8, "Add: Real-Real")
        r = a - b
        call check(mf_to_double(r) == 5.0_8, "Sub: Real-Real")
        r = a * b
        call check(mf_to_double(r) == 50.0_8, "Mul: Real-Real")
        r = a / b
        call check(mf_to_double(r) == 2.0_8, "Div: Real-Real")

        ! Real-DP
        r = a + d
        call check(mf_to_double(r) == 12.0_8, "Add: Real-DP")
        r = a - d
        call check(mf_to_double(r) == 8.0_8, "Sub: Real-DP")
        r = a * d
        call check(mf_to_double(r) == 20.0_8, "Mul: Real-DP")
        r = a / d
        call check(mf_to_double(r) == 5.0_8, "Div: Real-DP")

        ! DP-Real
        r = d + a
        call check(mf_to_double(r) == 12.0_8, "Add: DP-Real")
        r = d - a
        call check(mf_to_double(r) == -8.0_8, "Sub: DP-Real")
        r = d * a
        call check(mf_to_double(r) == 20.0_8, "Mul: DP-Real")
        r = d / a
        call check(mf_to_double(r) == 0.2_8, "Div: DP-Real")

        ! Real-Int
        r = a + i
        call check(mf_to_double(r) == 13.0_8, "Add: Real-Int")
        r = a - i
        call check(mf_to_double(r) == 7.0_8, "Sub: Real-Int")
        r = a * i
        call check(mf_to_double(r) == 30.0_8, "Mul: Real-Int")
        r = a / i
        call check(abs(mf_to_double(r) - 3.33333333333333_8) < 1e-14, "Div: Real-Int")

        ! Int-Real
        r = i + a
        call check(mf_to_double(r) == 13.0_8, "Add: Int-Real")
        r = i - a
        call check(mf_to_double(r) == -7.0_8, "Sub: Int-Real")
        ! r = i * a ! Missing Int-Real Mul in interface? Let me check.
        ! Actually I have int_mul_mf

        ! Unary and Power
        r = -a
        call check(mf_to_double(r) == -10.0_8, "Unary Minus")
        r = b**2
        call check(mf_to_double(r) == 25.0_8, "Power: Real**Int")
    end subroutine

    subroutine test_comparisons()
        type(float64x2) :: a, b
        real(8) :: d = 5.0_8

        a = float64x2(10.0_8)
        b = float64x2(5.0_8)

        call check(a > b, "GT: Real-Real")
        call check(b < a, "LT: Real-Real")
        call check(a == a, "EQ: Real-Real")
        call check(a /= b, "NE: Real-Real")
        call check(a >= b, "GE: Real-Real")
        call check(b <= a, "LE: Real-Real")

        call check(a > d, "GT: Real-DP")
        call check(b == d, "EQ: Real-DP")
        call check(d < a, "LT: DP-Real")
    end subroutine

    subroutine test_math_generics()
        type(float64x2) :: a, b, c, r
        
        a = float64x2(-1.0_8)
        call check(mf_to_double(abs(a)) == 1.0_8, "Math: ABS")
        
        a = float64x2(4.0_8)
        call check(mf_to_double(sqrt(a)) == 2.0_8, "Math: SQRT")

        a = float64x2(0.0_8)
        call check(mf_to_double(sin(a)) == 0.0_8, "Math: SIN")
        call check(mf_to_double(cos(a)) == 1.0_8, "Math: COS")

        a = float64x2(1.0_8)
        call check(abs(mf_to_double(exp(a)) - 2.71828182845904_8) < 1e-14, "Math: EXP")
        call check(mf_to_double(log(exp(a))) == 1.0_8, "Math: LOG")

        a = float64x2(100.0_8)
        call check(mf_to_double(log10(a)) == 2.0_8, "Math: LOG10")

        a = float64x2(1.0_8); b = float64x2(2.0_8); c = float64x2(3.0_8)
        call check(mf_to_double(min(a, b)) == 1.0_8, "Math: MIN(2)")
        call check(mf_to_double(min(a, b, c)) == 1.0_8, "Math: MIN(3)")
        call check(mf_to_double(max(a, b)) == 2.0_8, "Math: MAX(2)")
        call check(mf_to_double(max(a, b, c)) == 3.0_8, "Math: MAX(3)")

        call check(mf_to_double(sign(float64x2(5.0_8), float64x2(-1.0_8))) == -5.0_8, "Math: SIGN")
        call check(mf_to_double(mod(float64x2(7.0_8), float64x2(3.0_8))) == 1.0_8, "Math: MOD")
    end subroutine

    subroutine test_complex()
        type(float64x2) :: r, i
        type(complex128x2) :: z, w, res
        
        r = float64x2(3.0_8)
        i = float64x2(4.0_8)
        z = complex128x2(r, i)

        call check(mf_to_double(real(z)) == 3.0_8, "Complex: REAL part")
        call check(mf_to_double(aimag(z)) == 4.0_8, "Complex: AIMAG part")
        call check(mf_to_double(abs(z)) == 5.0_8, "Complex: ABS (magnitude)")

        w = conjg(z)
        call check(mf_to_double(aimag(w)) == -4.0_8, "Complex: CONJG")

        res = z + z
        call check(mf_to_double(real(res)) == 6.0_8 .and. mf_to_double(aimag(res)) == 8.0_8, "Complex: ADD")
    end subroutine

    subroutine test_io()
        type(float64x2) :: a, b
        character(len=100) :: buffer
        
        a = float64x2(1.23456789_8)
        
        ! Test write with DT and specific sub-format
        write(buffer, '(DT"F12.8")') a
        call check(adjustl(trim(buffer)) == "1.23456789", "Defined I/O: WRITE (DT)")

        ! Test read (usually uses default/list-directed logic in implementation)
        read(buffer, '(DT)') b
        call check(abs(mf_to_double(a) - mf_to_double(b)) < 1e-8, "Defined I/O: READ (DT)")

        ! Test list-directed
        print *, "List-directed output test (should show ~1.234): ", a
    end subroutine

end program multifloats_unittest
