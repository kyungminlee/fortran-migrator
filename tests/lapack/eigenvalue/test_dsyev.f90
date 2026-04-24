! Compares eigenvalue vector W (unique, ascending-sorted).
! Eigenvectors Z have a sign convention that's precision-dependent
! (a near-zero sign comparison inside LAPACK can flip), so we
! validate them via the invariant ||A*Z - Z*diag(W)|| / ||A|| which
! is immune to column-sign flips.
program test_dsyev
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsyev
    use ref_quad_lapack, only: dsyev
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    integer :: i, n, j, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: w_ref(:), w_got(:), work_ref(:)
    real(ep), allocatable :: AZ(:,:), ZD(:,:)
    real(ep) :: wopt(1), err_w, err_r, tol, anorm, rnorm
    character(len=48) :: label

    call report_init('dsyev', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 10001 + 71 * i)

        allocate(A_ref(n, n), A_got(n, n), w_ref(n), w_got(n))
        A_ref = A0;  A_got = A0

        call dsyev('V', 'U', n, A_ref, n, w_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work_ref(lwork))
        call dsyev('V', 'U', n, A_ref, n, w_ref, work_ref, lwork, info)

        call target_dsyev('V', 'U', n, A_got, n, w_got, info)

        err_w = max_rel_err_vec(w_got, w_ref)
        ! Invariant residual: ||A0 * Z_got - Z_got * diag(w_got)||.
        allocate(AZ(n, n), ZD(n, n))
        AZ = matmul(A0, A_got)
        do j = 1, n
            ZD(:, j) = A_got(:, j) * w_got(j)
        end do
        anorm = maxval(abs(A0))
        rnorm = maxval(abs(AZ - ZD))
        err_r = rnorm / max(anorm, tiny(1.0_ep))
        tol   = 16.0_ep * real(n, ep)**3 * target_eps

        write(label, '(a,i0,a)') 'n=', n, ',out=W'
        call report_case(trim(label), err_w, tol)
        write(label, '(a,i0,a)') 'n=', n, ',out=residual'
        call report_case(trim(label), err_r, tol)

        deallocate(A_ref, A_got, w_ref, w_got, work_ref, AZ, ZD)
    end do
    call report_finalize()
end program test_dsyev
