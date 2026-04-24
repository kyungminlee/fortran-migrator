program test_dtrmm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dtrmm
    use ref_quad_blas, only: dtrmm
    implicit none

    integer, parameter :: ms(*) = [4, 32, 100]
    integer, parameter :: ns(*) = [5, 40, 80]
    integer :: i, m, n
    real(ep), allocatable :: A(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: alpha, err, tol
    character(len=48) :: label

    call report_init('dtrmm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        ! Side='L': A is m-by-m triangular.
        call gen_matrix_quad(m, m, A,  seed = 911 + 23 * i)
        call gen_matrix_quad(m, n, B0, seed = 921 + 23 * i)
        alpha = real(0.7_ep, ep)
        allocate(B_ref(m, n), B_got(m, n))
        B_ref = B0
        B_got = B0
        call dtrmm('L', 'U', 'N', 'N', m, n, alpha, A, m, B_ref, m)
        call target_dtrmm('L', 'U', 'N', 'N', m, n, alpha, A, m, B_got, m)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * 2.0_ep * real(m, ep) * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(B_ref, B_got)
    end do
    call report_finalize()
end program test_dtrmm
