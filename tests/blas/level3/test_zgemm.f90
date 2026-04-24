program test_zgemm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex
    use target_blas,   only: target_name, target_eps, target_zgemm
    use ref_quad_blas, only: zgemm
    implicit none

    integer, parameter :: ms(*) = [4, 32, 80]
    integer, parameter :: ns(*) = [5, 40, 60]
    integer, parameter :: ks(*) = [6, 24, 50]
    integer :: i, m, n, k
    complex(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgemm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); k = ks(i)
        call gen_matrix_complex(m, k, A,  seed = 1001 + 29 * i)
        call gen_matrix_complex(k, n, B,  seed = 1011 + 29 * i)
        call gen_matrix_complex(m, n, C0, seed = 1021 + 29 * i)
        alpha = cmplx(0.7_ep, 0.2_ep, ep)
        beta  = cmplx(0.3_ep, -0.1_ep, ep)
        allocate(C_ref(m, n), C_got(m, n))
        C_ref = C0
        C_got = C0
        call zgemm('N', 'N', m, n, k, alpha, A, m, B, k, beta, C_ref, m)
        call target_zgemm('N', 'N', m, n, k, alpha, A, m, B, k, beta, C_got, m)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 32.0_ep * 8.0_ep * real(k, ep) * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(C_ref, C_got)
    end do
    call report_finalize()
end program test_zgemm
