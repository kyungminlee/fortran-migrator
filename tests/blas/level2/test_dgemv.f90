program test_dgemv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dgemv
    use ref_quad_blas, only: dgemv
    implicit none

    integer, parameter :: ms(*) = [5, 50, 200]
    integer, parameter :: ns(*) = [7, 60, 250]
    integer :: i, m, n
    real(ep), allocatable :: A(:,:), x(:), y0(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    character(len=32) :: label

    call report_init('dgemv', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A, seed = 401 + 17 * i)
        call gen_vector_quad(n,    x, seed = 411 + 17 * i)
        call gen_vector_quad(m,   y0, seed = 421 + 17 * i)
        alpha = real(0.6_ep, ep)
        beta  = real(0.3_ep, ep)
        allocate(y_ref(m), y_got(m))
        y_ref = y0
        y_got = y0
        call dgemv('N', m, n, alpha, A, m, x, 1, beta, y_ref, 1)
        call target_dgemv('N', m, n, alpha, A, m, x, 1, beta, y_got, 1)
        err = max_rel_err_vec(y_got, y_ref)
        tol = 16.0_ep * 2.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_dgemv
