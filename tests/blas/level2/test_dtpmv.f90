program test_dtpmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dtpmv
    use ref_quad_blas, only: dtpmv
    implicit none

    integer, parameter :: cases(*) = [10, 50, 200]
    integer :: i, n, ap_size
    real(ep), allocatable :: ap(:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('dtpmv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        ap_size = n * (n + 1) / 2
        call gen_vector_quad(ap_size, ap, seed = 621 + 17 * i)
        call gen_vector_quad(n,       x0, seed = 631 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0
        x_got = x0
        call dtpmv('U', 'N', 'N', n, ap, x_ref, 1)
        call target_dtpmv('U', 'N', 'N', n, ap, x_got, 1)
        err = max_rel_err_vec(x_got, x_ref)
        tol = 16.0_ep * 2.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_dtpmv
