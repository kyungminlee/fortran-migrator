program test_zgelq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgelq
    use ref_quad_lapack, only: zgelq
    implicit none

    integer, parameter :: ms(*) = [12, 20, 28]
    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, m, n, info, lwork, tsize, j, k
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    complex(ep), allocatable :: T_ref(:), T_got(:), L_ref(:,:), L_got(:,:)
    complex(ep) :: tquery(5), wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgelq', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 370031 + 47 * i)
        allocate(A_ref(m, n), A_got(m, n))
        A_ref = A0; A_got = A0
        call zgelq(m, n, A_ref, m, tquery, -1, wopt, -1, info)
        tsize = max(1, int(real(tquery(1), ep)))
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(T_ref(tsize), T_got(tsize), work(lwork))
        call zgelq(m, n, A_ref, m, T_ref, tsize, work, lwork, info)
        call target_zgelq(m, n, A_got, m, T_got, tsize, info)
        allocate(L_ref(m, m), L_got(m, m))
        L_ref = (0.0_ep, 0.0_ep); L_got = (0.0_ep, 0.0_ep)
        do k = 1, m
            do j = k, m
                L_ref(j, k) = cmplx(abs(A_ref(j, k)), 0.0_ep, ep)
                L_got(j, k) = cmplx(abs(A_got(j, k)), 0.0_ep, ep)
            end do
        end do
        err = max_rel_err_mat_z(L_got, L_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, T_ref, T_got, work, L_ref, L_got)
    end do
    call report_finalize()
end program test_zgelq
