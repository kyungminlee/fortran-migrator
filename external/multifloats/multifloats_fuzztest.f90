program multifloats_fuzztest
    use multifloats
    implicit none

    integer, parameter :: qp = 16
    integer, parameter :: dp = 8
    integer, parameter :: sp = 4
    integer, parameter :: num_iterations = 100000

    print *, "Starting Multifloats Fuzz Test (", num_iterations, " iterations)..."

    call fuzz_arithmetic_mf_mf()
    call fuzz_arithmetic_mixed()
    call fuzz_math_functions()
    call fuzz_complex_arithmetic()
    call fuzz_power()

    print *, "========================================================"
    print *, "Fuzz test complete."

contains

    ! ----------------------------------------------------------------
    ! Helpers
    ! ----------------------------------------------------------------

    function to_mf(v) result(res)
        real(qp), intent(in) :: v
        type(float64x2) :: res
        res%limbs(1) = real(v, dp)
        res%limbs(2) = real(v - real(res%limbs(1), qp), dp)
    end function

    function from_mf(mf) result(res)
        type(float64x2), intent(in) :: mf
        real(qp) :: res
        res = real(mf%limbs(1), qp) + real(mf%limbs(2), qp)
    end function

    function from_cx(cx) result(res)
        type(complex128x2), intent(in) :: cx
        complex(qp) :: res
        res = cmplx(from_mf(cx%re), from_mf(cx%im), kind=qp)
    end function

    subroutine update_stats(native, mf, max_e, avg_e)
        real(qp), intent(in) :: native
        type(float64x2), intent(in) :: mf
        real(qp), intent(inout) :: max_e, avg_e
        real(qp) :: err, val_mf

        val_mf = from_mf(mf)
        if (abs(native) > 0.0_qp) then
            err = abs(native - val_mf) / abs(native)
        else
            err = abs(native - val_mf)
        end if

        if (err > max_e) max_e = err
        avg_e = avg_e + err
    end subroutine

    subroutine update_stats_cx(native, cx, max_e, avg_e)
        complex(qp), intent(in) :: native
        type(complex128x2), intent(in) :: cx
        real(qp), intent(inout) :: max_e, avg_e
        real(qp) :: err
        complex(qp) :: val_cx

        val_cx = from_cx(cx)
        if (abs(native) > 0.0_qp) then
            err = abs(native - val_cx) / abs(native)
        else
            err = abs(native - val_cx)
        end if

        if (err > max_e) max_e = err
        avg_e = avg_e + err
    end subroutine

    subroutine print_stat(label, max_e, avg_e)
        character(*), intent(in) :: label
        real(qp), intent(in) :: max_e, avg_e
        print '(A12, E20.10, E20.10)', label, real(max_e, 8), real(avg_e, 8)
    end subroutine

    subroutine print_header(section)
        character(*), intent(in) :: section
        print *, "--------------------------------------------------------"
        print *, section
        print '(A12, A20, A20)', "Op", "Max Error", "Avg Error"
    end subroutine

    ! ----------------------------------------------------------------
    ! Fuzz: mf OP mf
    ! ----------------------------------------------------------------

    subroutine fuzz_arithmetic_mf_mf()
        real(qp) :: max_add, avg_add, max_sub, avg_sub
        real(qp) :: max_mul, avg_mul, max_div, avg_div
        real(qp) :: v1, v2
        type(float64x2) :: mf1, mf2
        real(8) :: rnd(2)
        integer :: i

        max_add = 0; avg_add = 0; max_sub = 0; avg_sub = 0
        max_mul = 0; avg_mul = 0; max_div = 0; avg_div = 0

        do i = 1, num_iterations
            call random_number(rnd)
            v1 = real(rnd(1), qp) * 100.0_qp + 0.1_qp
            call random_number(rnd)
            v2 = real(rnd(2), qp) * 100.0_qp + 0.1_qp

            mf1 = to_mf(v1)
            mf2 = to_mf(v2)

            call update_stats(v1 + v2, mf1 + mf2, max_add, avg_add)
            call update_stats(v1 - v2, mf1 - mf2, max_sub, avg_sub)
            call update_stats(v1 * v2, mf1 * mf2, max_mul, avg_mul)
            call update_stats(v1 / v2, mf1 / mf2, max_div, avg_div)
        end do

        call print_header("mf OP mf")
        call print_stat("ADD", max_add, avg_add / num_iterations)
        call print_stat("SUB", max_sub, avg_sub / num_iterations)
        call print_stat("MUL", max_mul, avg_mul / num_iterations)
        call print_stat("DIV", max_div, avg_div / num_iterations)
    end subroutine

    ! ----------------------------------------------------------------
    ! Fuzz: mixed-type arithmetic (mf OP dp, dp OP mf, mf OP int, etc.)
    ! ----------------------------------------------------------------

    subroutine fuzz_arithmetic_mixed()
        real(qp) :: max_mf_dp, avg_mf_dp, max_dp_mf, avg_dp_mf
        real(qp) :: max_mf_sp, avg_mf_sp, max_sp_mf, avg_sp_mf
        real(qp) :: max_mf_int, avg_mf_int, max_int_mf, avg_int_mf
        real(qp) :: v1
        real(8) :: d
        real(4) :: s
        integer :: iv
        type(float64x2) :: mf1
        real(8) :: rnd(2)
        integer :: i

        max_mf_dp = 0; avg_mf_dp = 0; max_dp_mf = 0; avg_dp_mf = 0
        max_mf_sp = 0; avg_mf_sp = 0; max_sp_mf = 0; avg_sp_mf = 0
        max_mf_int = 0; avg_mf_int = 0; max_int_mf = 0; avg_int_mf = 0

        do i = 1, num_iterations
            call random_number(rnd)
            v1 = real(rnd(1), qp) * 100.0_qp + 0.1_qp
            d = rnd(2) * 100.0_8 + 0.1_8
            s = real(rnd(2) * 100.0, sp) + 0.1
            iv = int(rnd(2) * 100.0) + 1

            mf1 = to_mf(v1)

            ! mf * dp
            call update_stats(v1 * real(d, qp), mf1 * d, max_mf_dp, avg_mf_dp)
            ! dp * mf
            call update_stats(real(d, qp) * v1, d * mf1, max_dp_mf, avg_dp_mf)
            ! mf * sp
            call update_stats(v1 * real(s, qp), mf1 * s, max_mf_sp, avg_mf_sp)
            ! sp * mf
            call update_stats(real(s, qp) * v1, s * mf1, max_sp_mf, avg_sp_mf)
            ! mf + int
            call update_stats(v1 + real(iv, qp), mf1 + iv, max_mf_int, avg_mf_int)
            ! int + mf
            call update_stats(real(iv, qp) + v1, iv + mf1, max_int_mf, avg_int_mf)
        end do

        call print_header("Mixed-type arithmetic (multiply/add)")
        call print_stat("mf * dp", max_mf_dp, avg_mf_dp / num_iterations)
        call print_stat("dp * mf", max_dp_mf, avg_dp_mf / num_iterations)
        call print_stat("mf * sp", max_mf_sp, avg_mf_sp / num_iterations)
        call print_stat("sp * mf", max_sp_mf, avg_sp_mf / num_iterations)
        call print_stat("mf + int", max_mf_int, avg_mf_int / num_iterations)
        call print_stat("int + mf", max_int_mf, avg_int_mf / num_iterations)
    end subroutine

    ! ----------------------------------------------------------------
    ! Fuzz: math functions
    ! ----------------------------------------------------------------

    subroutine fuzz_math_functions()
        real(qp) :: max_sqrt, avg_sqrt, max_exp, avg_exp
        real(qp) :: max_sin, avg_sin, max_cos, avg_cos
        real(qp) :: max_atan, avg_atan, max_atan2, avg_atan2
        real(qp) :: max_asin, avg_asin, max_log, avg_log
        real(qp) :: v1, v2
        type(float64x2) :: mf1, mf2
        real(8) :: rnd(2)
        integer :: i

        max_sqrt = 0; avg_sqrt = 0; max_exp = 0; avg_exp = 0
        max_sin = 0; avg_sin = 0; max_cos = 0; avg_cos = 0
        max_atan = 0; avg_atan = 0; max_atan2 = 0; avg_atan2 = 0
        max_asin = 0; avg_asin = 0; max_log = 0; avg_log = 0

        do i = 1, num_iterations
            call random_number(rnd)
            v1 = real(rnd(1), qp) * 100.0_qp + 0.1_qp
            v2 = real(rnd(2), qp) * 100.0_qp + 0.1_qp

            mf1 = to_mf(v1)
            mf2 = to_mf(v2)

            call update_stats(sqrt(v1), sqrt(mf1), max_sqrt, avg_sqrt)
            call update_stats(log(v1), log(mf1), max_log, avg_log)
            call update_stats(exp(real(rnd(1), qp)), exp(to_mf(real(rnd(1), qp))), max_exp, avg_exp)
            call update_stats(sin(v1), sin(mf1), max_sin, avg_sin)
            call update_stats(cos(v1), cos(mf1), max_cos, avg_cos)
            call update_stats(atan(v1), atan(mf1), max_atan, avg_atan)
            call update_stats(atan2(v1, v2), atan2(mf1, mf2), max_atan2, avg_atan2)

            ! asin needs input in [-1, 1]
            block
                real(qp) :: asin_v
                asin_v = real(rnd(1), qp) * 2.0_qp - 1.0_qp  ! [-1, 1]
                call update_stats(asin(asin_v), asin(to_mf(asin_v)), max_asin, avg_asin)
            end block
        end do

        call print_header("Math functions")
        call print_stat("SQRT", max_sqrt, avg_sqrt / num_iterations)
        call print_stat("LOG", max_log, avg_log / num_iterations)
        call print_stat("EXP", max_exp, avg_exp / num_iterations)
        call print_stat("SIN", max_sin, avg_sin / num_iterations)
        call print_stat("COS", max_cos, avg_cos / num_iterations)
        call print_stat("ATAN", max_atan, avg_atan / num_iterations)
        call print_stat("ATAN2", max_atan2, avg_atan2 / num_iterations)
        call print_stat("ASIN", max_asin, avg_asin / num_iterations)
    end subroutine

    ! ----------------------------------------------------------------
    ! Fuzz: complex arithmetic
    ! ----------------------------------------------------------------

    subroutine fuzz_complex_arithmetic()
        real(qp) :: max_add, avg_add, max_mul, avg_mul, max_div, avg_div
        real(qp) :: max_abs, avg_abs
        complex(qp) :: c1, c2
        type(complex128x2) :: z1, z2
        real(8) :: rnd(4)
        integer :: i

        max_add = 0; avg_add = 0; max_mul = 0; avg_mul = 0
        max_div = 0; avg_div = 0; max_abs = 0; avg_abs = 0

        do i = 1, num_iterations
            call random_number(rnd)
            c1 = cmplx(rnd(1)*100.0_8 + 0.1, rnd(2)*100.0_8 + 0.1, kind=qp)
            c2 = cmplx(rnd(3)*100.0_8 + 0.1, rnd(4)*100.0_8 + 0.1, kind=qp)

            z1%re = to_mf(real(c1, qp))
            z1%im = to_mf(aimag(c1))
            z2%re = to_mf(real(c2, qp))
            z2%im = to_mf(aimag(c2))

            call update_stats_cx(c1 + c2, z1 + z2, max_add, avg_add)
            call update_stats_cx(c1 * c2, z1 * z2, max_mul, avg_mul)
            call update_stats_cx(c1 / c2, z1 / z2, max_div, avg_div)

            ! abs(cx) → mf (compare as real)
            call update_stats(abs(c1), abs(z1), max_abs, avg_abs)
        end do

        call print_header("Complex arithmetic")
        call print_stat("CX ADD", max_add, avg_add / num_iterations)
        call print_stat("CX MUL", max_mul, avg_mul / num_iterations)
        call print_stat("CX DIV", max_div, avg_div / num_iterations)
        call print_stat("CX ABS", max_abs, avg_abs / num_iterations)
    end subroutine

    ! ----------------------------------------------------------------
    ! Fuzz: power
    ! ----------------------------------------------------------------

    subroutine fuzz_power()
        real(qp) :: max_pow_int, avg_pow_int, max_pow_mf, avg_pow_mf
        real(qp) :: v1
        type(float64x2) :: mf1
        real(8) :: rnd(2)
        integer :: i, ipow

        max_pow_int = 0; avg_pow_int = 0
        max_pow_mf = 0; avg_pow_mf = 0

        do i = 1, num_iterations
            call random_number(rnd)
            v1 = real(rnd(1), qp) * 10.0_qp + 0.1_qp  ! Keep base modest
            ipow = int(rnd(2) * 8.0)  ! 0..7

            mf1 = to_mf(v1)

            ! mf ** int
            call update_stats(v1 ** ipow, mf1 ** ipow, max_pow_int, avg_pow_int)

            ! mf ** mf (use small exponent to avoid overflow)
            block
                real(qp) :: exp_v
                exp_v = real(rnd(2), qp) * 3.0_qp
                call update_stats(v1 ** exp_v, mf1 ** to_mf(exp_v), max_pow_mf, avg_pow_mf)
            end block
        end do

        call print_header("Power")
        call print_stat("mf ** int", max_pow_int, avg_pow_int / num_iterations)
        call print_stat("mf ** mf", max_pow_mf, avg_pow_mf / num_iterations)
    end subroutine

end program multifloats_fuzztest
