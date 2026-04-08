program multifloats_fuzztest
    use multifloats
    implicit none

    integer, parameter :: qp = 16
    integer, parameter :: dp = 8
    integer, parameter :: num_iterations = 100000
    
    real(qp) :: max_err_add = 0.0_qp, avg_err_add = 0.0_qp
    real(qp) :: max_err_sub = 0.0_qp, avg_err_sub = 0.0_qp
    real(qp) :: max_err_mul = 0.0_qp, avg_err_mul = 0.0_qp
    real(qp) :: max_err_div = 0.0_qp, avg_err_div = 0.0_qp
    real(qp) :: max_err_sqrt = 0.0_qp, avg_err_sqrt = 0.0_qp
    real(qp) :: max_err_exp = 0.0_qp, avg_err_exp = 0.0_qp
    real(qp) :: max_err_sin = 0.0_qp, avg_err_sin = 0.0_qp

    integer :: i
    real(qp) :: v1, v2, res_native, res_mf_qp
    type(float64x2) :: mf1, mf2, mf_res
    real(8) :: rnd(2)

    print *, "Starting Multifloats Fuzz Test (", num_iterations, " iterations)..."

    do i = 1, num_iterations
        ! Generate random values in [0, 1] then scale
        call random_number(rnd)
        v1 = real(rnd(1), qp) * 100.0_qp + 0.1_qp  ! Avoid zero
        call random_number(rnd)
        v2 = real(rnd(2), qp) * 100.0_qp + 0.1_qp

        mf1 = to_mf(v1)
        mf2 = to_mf(v2)

        ! Addition
        res_native = v1 + v2
        mf_res = mf1 + mf2
        call update_stats(res_native, mf_res, max_err_add, avg_err_add)

        ! Subtraction
        res_native = v1 - v2
        mf_res = mf1 - mf2
        call update_stats(res_native, mf_res, max_err_sub, avg_err_sub)

        ! Multiplication
        res_native = v1 * v2
        mf_res = mf1 * mf2
        call update_stats(res_native, mf_res, max_err_mul, avg_err_mul)

        ! Division
        res_native = v1 / v2
        mf_res = mf1 / mf2
        call update_stats(res_native, mf_res, max_err_div, avg_err_div)

        ! Sqrt
        res_native = sqrt(v1)
        mf_res = sqrt(mf1)
        call update_stats(res_native, mf_res, max_err_sqrt, avg_err_sqrt)

        ! Exp
        res_native = exp(real(rnd(1), qp)) ! Keep in safe range
        mf_res = exp(to_mf(real(rnd(1), qp)))
        call update_stats(res_native, mf_res, max_err_exp, avg_err_exp)

        ! Sin
        res_native = sin(v1)
        mf_res = sin(mf1)
        call update_stats(res_native, mf_res, max_err_sin, avg_err_sin)
    end do

    avg_err_add = avg_err_add / num_iterations
    avg_err_sub = avg_err_sub / num_iterations
    avg_err_mul = avg_err_mul / num_iterations
    avg_err_div = avg_err_div / num_iterations
    avg_err_sqrt = avg_err_sqrt / num_iterations
    avg_err_exp = avg_err_exp / num_iterations
    avg_err_sin = avg_err_sin / num_iterations

    print *, "--------------------------------------------------------"
    print *, "Relative Error Statistics (compared to REAL(16)):"
    print '(A10, A20, A20)', "Op", "Max Error", "Avg Error"
    call print_stat("ADD", max_err_add, avg_err_add)
    call print_stat("SUB", max_err_sub, avg_err_sub)
    call print_stat("MUL", max_err_mul, avg_err_mul)
    call print_stat("DIV", max_err_div, avg_err_div)
    call print_stat("SQRT", max_err_sqrt, avg_err_sqrt)
    call print_stat("EXP", max_err_exp, avg_err_exp)
    call print_stat("SIN", max_err_sin, avg_err_sin)
    print *, "--------------------------------------------------------"

contains

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

    subroutine print_stat(label, max_e, avg_e)
        character(*), intent(in) :: label
        real(qp), intent(in) :: max_e, avg_e
        print '(A10, E20.10, E20.10)', label, real(max_e, 8), real(avg_e, 8)
    end subroutine

end program multifloats_fuzztest
