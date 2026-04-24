program test_dtbmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dtbmv
    use ref_quad_blas, only: dtbmv
    implicit none

    integer, parameter :: cases(*) = [20, 100]
    integer :: i, n, k, lda
    real(ep), allocatable :: A(:,:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('dtbmv', target_name)
    do i = 1, size(cases)
        n   = cases(i)
        k   = 3
        lda = k + 1
        call gen_matrix_quad(lda, n, A,  seed = 601 + 17 * i)
        call gen_vector_quad(n,      x0, seed = 611 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0
        x_got = x0
        call dtbmv('U', 'N', 'N', n, k, A, lda, x_ref, 1)
        call target_dtbmv('U', 'N', 'N', n, k, A, lda, x_got, 1)
        err = max_rel_err_vec(x_got, x_ref)
        tol = 16.0_ep * 2.0_ep * real(lda, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_dtbmv
