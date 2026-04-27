program test_ztzpad
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat_z
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_ztzpad
    use test_data,            only: gen_matrix_complex
    use target_ptzblas,       only: target_name, target_eps, target_ztzpad
    implicit none

    integer, parameter :: m = 16, n = 16
    character(len=1), parameter :: uplo_set(*) = ['L', 'U', 'L', 'U']
    character(len=1), parameter :: herm_set(*) = ['N', 'N', 'Z', 'Z']
    integer :: i
    character :: uplo, herm
    complex(ep), allocatable :: A_got(:,:), A_ref(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztzpad', target_name, 0)
    do i = 1, size(uplo_set)
        uplo = uplo_set(i); herm = herm_set(i)
        call gen_matrix_complex(m, n, A_got, seed = 3300 + 7 * i)
        A_ref = A_got
        alpha = cmplx(0.25_ep + 0.1_ep * i, 0.3_ep, ep)
        beta  = cmplx(-1.0_ep, 0.5_ep, ep)
        call target_ztzpad(uplo, herm, m, n, 0, alpha, beta, A_got, m)
        call ref_ztzpad(uplo, herm, m, n, 0, alpha, beta, A_ref)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 4.0_ep * target_eps
        write(label, '(a,a,a,a)') 'uplo=', uplo, ',herm=', herm
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_ztzpad
