program test_dgelqt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgelqt
    use ref_quad_lapack, only: dgelqt
    implicit none

    integer, parameter :: ms(*) = [12, 20, 28]
    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: mb    = 4
    integer :: i, m, n, info, j, k
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), T_ref(:,:), T_got(:,:), work(:)
    real(ep), allocatable :: L_ref(:,:), L_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgelqt', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 380021 + 47 * i)
        allocate(A_ref(m, n), A_got(m, n), T_ref(mb, min(m, n)), T_got(mb, min(m, n)), work(mb*m))
        A_ref = A0; A_got = A0
        call dgelqt(m, n, mb, A_ref, m, T_ref, mb, work, info)
        call target_dgelqt(m, n, mb, A_got, m, T_got, mb, info)
        allocate(L_ref(m, m), L_got(m, m))
        L_ref = 0.0_ep; L_got = 0.0_ep
        do k = 1, m
            do j = k, m
                L_ref(j, k) = abs(A_ref(j, k))
                L_got(j, k) = abs(A_got(j, k))
            end do
        end do
        err = max_rel_err_mat(L_got, L_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, T_ref, T_got, work, L_ref, L_got)
    end do
    call report_finalize()
end program test_dgelqt
