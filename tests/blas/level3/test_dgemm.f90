program test_dgemm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dgemm
    use ref_quad_blas, only: dgemm
    implicit none

    integer, parameter :: ms(*) = [4, 32, 128, 16, 16]
    integer, parameter :: ns(*) = [5, 40, 100, 12, 12]
    integer, parameter :: ks(*) = [6, 24,  80,  8,  8]
    ! Per-shape (alpha, beta) — the last two cases exercise BLAS's
    ! quick-return paths: ``alpha=0`` (skip the matmul; just scale C
    ! by beta) and ``beta=0`` (overwrite C without reading initial
    ! contents — a path that bypasses entirely separate code in
    ! every blocked GEMM kernel).
    real(ep), parameter :: alphas(*) = [0.7_ep, 0.7_ep, 0.7_ep, 0.0_ep, 0.5_ep]
    real(ep), parameter :: betas(*)  = [0.3_ep, 0.3_ep, 0.3_ep, 0.4_ep, 0.0_ep]
    integer :: i, m, n, k
    real(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call report_init('dgemm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); k = ks(i)
        call gen_matrix_quad(m, k, A,  seed = 801 + 23 * i)
        call gen_matrix_quad(k, n, B,  seed = 811 + 23 * i)
        call gen_matrix_quad(m, n, C0, seed = 821 + 23 * i)
        alpha = alphas(i); beta = betas(i)
        allocate(C_ref(m, n), C_got(m, n))
        C_ref = C0
        C_got = C0
        call dgemm('N', 'N', m, n, k, alpha, A, m, B, k, beta, C_ref, m)
        call target_dgemm('N', 'N', m, n, k, alpha, A, m, B, k, beta, C_got, m)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * 2.0_ep * real(k, ep) * target_eps
        write(label, '(a,f3.1,a,f3.1,a,i0,a,i0,a,i0)') &
            'a=', alpha, ',b=', beta, ',m=', m, ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, B, C0, C_ref, C_got)
    end do
    call report_finalize()
end program test_dgemm
