program test_ztrsm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex
    use target_blas,   only: target_name, target_eps, target_ztrsm
    use ref_quad_blas, only: ztrsm
    implicit none

    integer, parameter :: ms(*) = [4, 32, 80]
    integer, parameter :: ns(*) = [5, 40, 60]
    integer :: i, m, n, j
    complex(ep), allocatable :: A(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztrsm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, m, A,  seed = 1081 + 29 * i)
        ! Beef up diagonal so the triangle is well-conditioned.
        do j = 1, m
            A(j, j) = A(j, j) + cmplx(real(m + 2, ep), 0.0_ep, ep)
        end do
        call gen_matrix_complex(m, n, B0, seed = 1091 + 29 * i)
        alpha = cmplx(0.7_ep, 0.2_ep, ep)
        allocate(B_ref(m, n), B_got(m, n))
        B_ref = B0
        B_got = B0
        call ztrsm('L', 'U', 'N', 'N', m, n, alpha, A, m, B_ref, m)
        call target_ztrsm('L', 'U', 'N', 'N', m, n, alpha, A, m, B_got, m)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 32.0_ep * 8.0_ep * real(m, ep) * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(B_ref, B_got)
    end do
    call report_finalize()
end program test_ztrsm
