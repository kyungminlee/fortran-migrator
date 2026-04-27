! dhsein: inverse iteration eigenvectors of an upper Hessenberg matrix
! given approximate eigenvalues. SIDE='R', EIGSRC='N', INITV='N'.
program test_dhsein
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dhsein
    use ref_quad_lapack, only: dgehrd, dhseqr, dhsein
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, m_r, m_g
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: H(:,:), WR(:), WI(:), WR_r(:), Z(:,:)
    real(ep), allocatable :: VL(:,:), VR_r(:,:), VR_g(:,:)
    integer,  allocatable :: ifaill(:), ifailr(:)
    logical,  allocatable :: select(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dhsein', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 112001 + 47 * i)
        allocate(tau(n - 1))
        call dgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(H(n,n), WR(n), WI(n), WR_r(n), Z(1,1))
        H = A
        call dhseqr('E', 'N', n, 1, n, H, n, WR, WI, Z, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dhseqr('E', 'N', n, 1, n, H, n, WR, WI, Z, 1, work, lwork, info)
        deallocate(work)
        ! Restore H to its Hessenberg form (overwritten by dhseqr).
        H = A
        allocate(VL(1, n), VR_r(n, n), VR_g(n, n))
        allocate(ifaill(n), ifailr(n), select(n), work(n*(n + 2)))
        select = .true.
        VR_r = 0.0_ep; VR_g = 0.0_ep
        WR_r = WR
        call dhsein('R', 'N', 'N', select, n, H, n, WR_r, WI, &
                    VL, 1, VR_r, n, n, m_r, work, ifaill, ifailr, info)
        deallocate(work)
        WR_r = WR
        select = .true.
        call target_dhsein('R', 'N', 'N', select, n, H, n, WR_r, WI, &
                           VL, 1, VR_g, n, n, m_g, ifaill, ifailr, info)
        err = max_rel_err_mat(VR_g, VR_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, H, WR, WI, WR_r, Z, VL, VR_r, VR_g, &
                   ifaill, ifailr, select)
    end do
    call report_finalize()
end program test_dhsein
