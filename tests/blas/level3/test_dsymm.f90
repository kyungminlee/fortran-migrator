program test_dsymm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dsymm
    use ref_quad_blas, only: dsymm
    implicit none

    integer, parameter :: ms(*) = [4, 32, 100]
    integer, parameter :: ns(*) = [5, 40, 80]
    integer :: i, m, n
    real(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call report_init('dsymm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        ! Side='L': A is m-by-m symmetric.
        call gen_matrix_quad(m, m, A,  seed = 831 + 23 * i)
        call gen_matrix_quad(m, n, B,  seed = 841 + 23 * i)
        call gen_matrix_quad(m, n, C0, seed = 851 + 23 * i)
        alpha = real(0.6_ep, ep)
        beta  = real(0.4_ep, ep)
        allocate(C_ref(m, n), C_got(m, n))
        C_ref = C0
        C_got = C0
        call dsymm('L', 'U', m, n, alpha, A, m, B, m, beta, C_ref, m)
        call target_dsymm('L', 'U', m, n, alpha, A, m, B, m, beta, C_got, m)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * 2.0_ep * real(m, ep) * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(C_ref, C_got)
    end do
    call report_finalize()
end program test_dsymm
