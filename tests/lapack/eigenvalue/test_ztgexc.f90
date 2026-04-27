! ztgexc: complex generalized Schur exchange. Smoke test.
! KNOWN FAILING (Phase L4): SIGSEGV. Synthetic B=I setup may
! collide with ztgexc's input-validation. See tests/lapack/TODO.md.
program test_ztgexc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztgexc
    use ref_quad_lapack, only: zgees, ztgexc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, sdim, j
    complex(ep), allocatable :: A(:,:), B(:,:), AR(:,:), BR(:,:), Q(:,:), Z(:,:)
    complex(ep), allocatable :: AG(:,:), BG(:,:), QG(:,:), ZG(:,:), W(:), VS(:,:), work(:)
    real(ep),    allocatable :: rwork(:), m_r(:), m_g(:)
    logical,     allocatable :: bwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztgexc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 141001 + 47 * i)
        allocate(B(n,n)); B = (0.0_ep, 0.0_ep)
        do j = 1, n; B(j,j) = (1.0_ep, 0.0_ep); end do  ! identity B (problem reduces to standard)
        allocate(AR(n,n), BR(n,n), Q(n,n), Z(n,n), AG(n,n), BG(n,n), QG(n,n), ZG(n,n))
        allocate(W(n), VS(n,n), rwork(n), bwork(n), m_r(n), m_g(n))
        AR = A
        call zgees('V', 'N', sel_all_c, n, AR, n, sdim, W, VS, n, &
                   wopt, -1, rwork, bwork, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgees('V', 'N', sel_all_c, n, AR, n, sdim, W, VS, n, &
                   work, lwork, rwork, bwork, info)
        deallocate(work)
        BR = B
        Q = VS; Z = VS
        AG = AR; BG = BR; QG = Q; ZG = Z
        call ztgexc(.true., .true., n, AR, n, BR, n, Q, n, Z, n, n, 1, info)
        call target_ztgexc(.true., .true., n, AG, n, BG, n, QG, n, ZG, n, n, 1, info)
        do j = 1, n
            m_r(j) = abs(AR(j,j))
            m_g(j) = abs(AG(j,j))
        end do
        call sort_desc(m_r, n); call sort_desc(m_g, n)
        err = max_rel_err_vec(m_g, m_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, AR, BR, Q, Z, AG, BG, QG, ZG, W, VS, rwork, bwork, m_r, m_g)
    end do
    call report_finalize()
contains
    logical function sel_all_c(z)
        complex(ep), intent(in) :: z
        sel_all_c = .true.
        if (z == z) continue
    end function
    subroutine sort_desc(x, m)
        real(ep), intent(inout) :: x(:)
        integer,  intent(in)    :: m
        integer :: ii, jj
        real(ep) :: tt
        do ii = 1, m - 1
            do jj = ii + 1, m
                if (x(ii) < x(jj)) then
                    tt = x(ii); x(ii) = x(jj); x(jj) = tt
                end if
            end do
        end do
    end subroutine
end program test_ztgexc
