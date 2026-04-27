program test_dvasum
    use prec_kinds,            only: ep
    use compare,               only: rel_err_scalar
    use ptzblas_prec_report,   only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_blas, only: dasum
    use test_data,             only: gen_vector_quad
    use target_ptzblas,        only: target_name, target_eps, target_dvasum
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer :: i, n
    real(ep), allocatable :: x(:)
    real(ep) :: ref, got, err, tol
    character(len=32) :: label

    call report_init('dvasum', target_name, 0)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x, seed = 601 + 7 * i)
        got = 0.0_ep
        call target_dvasum(n, got, x, 1)
        ref = dasum(n, x, 1)
        err = rel_err_scalar(got, ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x)
    end do
    call report_finalize()
end program test_dvasum
