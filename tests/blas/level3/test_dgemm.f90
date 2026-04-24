program test_dgemm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dgemm
    use ref_quad_blas, only: dgemm
    implicit none

    integer, parameter :: ms(*) = [4, 32, 128]
    integer, parameter :: ns(*) = [5, 40, 100]
    integer, parameter :: ks(*) = [6, 24, 80]
    integer :: i, m, n, k
    real(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call report_init('dgemm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); k = ks(i)
        call gen_matrix_quad(m, k, A,  seed = 801 + 23 * i)
        call gen_matrix_quad(k, n, B,  seed = 811 + 23 * i)
        call gen_matrix_quad(m, n, C0, seed = 821 + 23 * i)
        alpha = real(0.7_ep, ep)
        beta  = real(0.3_ep, ep)
        allocate(C_ref(m, n), C_got(m, n))
        C_ref = C0
        C_got = C0
        call dgemm('N', 'N', m, n, k, alpha, A, m, B, k, beta, C_ref, m)
        call target_dgemm('N', 'N', m, n, k, alpha, A, m, B, k, beta, C_got, m)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * 2.0_ep * real(k, ep) * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(C_ref, C_got)
    end do
    call report_finalize()
end program test_dgemm
