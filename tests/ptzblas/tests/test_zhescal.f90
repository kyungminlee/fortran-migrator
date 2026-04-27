program test_zhescal
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat_z
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_zhescal
    use test_data,            only: gen_matrix_complex
    use target_ptzblas,       only: target_name, target_eps, target_zhescal
    implicit none

    integer, parameter :: m = 24, n = 24
    character(len=1), parameter :: uplo_set(*) = ['L', 'U']
    integer :: i
    character :: uplo
    complex(ep), allocatable :: A_got(:,:), A_ref(:,:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('zhescal', target_name, 0)
    do i = 1, size(uplo_set)
        uplo = uplo_set(i)
        call gen_matrix_complex(m, n, A_got, seed = 2700 + 5 * i)
        A_ref = A_got
        alpha = 1.4_ep + 0.1_ep * real(i, ep)
        call target_zhescal(uplo, m, n, 0, alpha, A_got, m)
        call ref_zhescal(uplo, m, n, 0, alpha, A_ref)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 8.0_ep * target_eps
        write(label, '(a,a)') 'uplo=', uplo
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_zhescal
