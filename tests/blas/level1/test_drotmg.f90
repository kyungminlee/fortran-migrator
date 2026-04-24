program test_drotmg
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar
    use target_blas,   only: target_name, target_eps, target_drotmg
    use ref_quad_blas, only: drotmg
    implicit none

    integer, parameter :: ncases = 4
    integer :: i
    real(ep) :: d1_0, d2_0, x1_0, y1_0
    real(ep) :: d1_ref, d2_ref, x1_ref, param_ref(5)
    real(ep) :: d1_got, d2_got, x1_got, param_got(5)
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('drotmg', target_name)
    do i = 1, ncases
        d1_0 = real(i, ep) * 1.5_ep
        d2_0 = real(ncases - i + 1, ep) * 0.8_ep
        x1_0 = real(i, ep) * 0.5_ep
        y1_0 = real(i, ep) * 0.3_ep
        d1_ref = d1_0; d2_ref = d2_0; x1_ref = x1_0
        d1_got = d1_0; d2_got = d2_0; x1_got = x1_0
        call drotmg(d1_ref, d2_ref, x1_ref, y1_0, param_ref)
        call target_drotmg(d1_got, d2_got, x1_got, y1_0, param_got)
        ! Compare d1, d2, x1, and the FLAG (param(1)). Other param
        ! slots are flag-conditional (BLAS spec: p(2)/p(5) only used
        ! when flag=1; p(3)/p(4) only when flag=0; etc.) — comparing
        ! the unused slots is a false positive.
        err = max(rel_err_scalar(d1_got, d1_ref), &
                  rel_err_scalar(d2_got, d2_ref), &
                  rel_err_scalar(x1_got, x1_ref), &
                  rel_err_scalar(param_got(1), param_ref(1)))
        call compare_meaningful_params(param_ref, param_got, err)
        tol = 32.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()

contains

    ! Compare only the param slots that are meaningful for the
    ! returned flag value. Updates `err` with the max relative error
    ! across the meaningful slots.
    subroutine compare_meaningful_params(p_ref, p_got, err)
        real(ep), intent(in)    :: p_ref(5), p_got(5)
        real(ep), intent(inout) :: err
        integer :: flag

        flag = nint(p_ref(1))
        select case (flag)
        case (-1)
            ! All four off-diagonal/diagonal entries are meaningful
            err = max(err, rel_err_scalar(p_got(2), p_ref(2)), &
                          rel_err_scalar(p_got(3), p_ref(3)), &
                          rel_err_scalar(p_got(4), p_ref(4)), &
                          rel_err_scalar(p_got(5), p_ref(5)))
        case (0)
            ! H = [1, h12; h21, 1] — only p(3)=h12, p(4)=h21 used
            err = max(err, rel_err_scalar(p_got(3), p_ref(3)), &
                          rel_err_scalar(p_got(4), p_ref(4)))
        case (1)
            ! H = [h11, 1; -1, h22] — only p(2)=h11, p(5)=h22 used
            err = max(err, rel_err_scalar(p_got(2), p_ref(2)), &
                          rel_err_scalar(p_got(5), p_ref(5)))
        case (-2)
            ! H = identity — no off-diagonal params used
        end select
    end subroutine compare_meaningful_params

end program test_drotmg
