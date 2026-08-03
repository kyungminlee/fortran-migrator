program test_pzgemm
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zgemm
    use pblas_grid,    only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib, only: gen_distrib_matrix_z, gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pzgemm
    implicit none

    ! Cases:
    !   1-3: TRANSA=N, TRANSB=N — main path
    !   4:   alpha=0 quick-return (parity with pdgemm)
    !   5:   beta=0 quick-return
    !   6:   TRANSA=T, TRANSB=N
    !   7:   TRANSA=N, TRANSB=T
    !   8:   TRANSA=C, TRANSB=N — conjugate-transpose path
    !   9:   TRANSA=N, TRANSB=C
    integer, parameter :: ms(*) = [24, 64, 96, 32, 32, 48, 48, 48, 48]
    integer, parameter :: ns(*) = [32, 48, 80, 32, 32, 40, 40, 40, 40]
    integer, parameter :: ks(*) = [16, 40, 72, 24, 24, 32, 32, 32, 32]
    character(len=1), parameter :: transas(*) = ['N','N','N','N','N','T','N','C','N']
    character(len=1), parameter :: transbs(*) = ['N','N','N','N','N','N','T','N','C']
    real(ep), parameter :: alpha_re(*) = [0.7_ep, 0.7_ep, 0.7_ep, 0.0_ep, &
                                           0.5_ep, 0.6_ep, 0.6_ep, 0.6_ep, 0.6_ep]
    real(ep), parameter :: alpha_im(*) = [0.2_ep, 0.2_ep, 0.2_ep, 0.0_ep, &
                                           0.1_ep, 0.3_ep, 0.3_ep, 0.3_ep, 0.3_ep]
    real(ep), parameter :: beta_re(*)  = [0.3_ep, 0.3_ep, 0.3_ep, 0.4_ep, &
                                           0.0_ep, 0.2_ep, 0.2_ep, 0.2_ep, 0.2_ep]
    real(ep), parameter :: beta_im(*)  = [-0.1_ep, -0.1_ep, -0.1_ep, 0.1_ep, &
                                            0.0_ep, -0.05_ep, -0.05_ep, -0.05_ep, -0.05_ep]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, k
    integer :: ar, ac, br, bc
    integer :: desca(9), descb(9), descc(9)
    character(len=1) :: transa, transb
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), B_glob(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call grid_init()
    call report_init('pzgemm', target_name, my_rank)

    do i = 1, size(ms)
        m = ms(i); n = ns(i); k = ks(i)
        alpha = cmplx(alpha_re(i), alpha_im(i), ep)
        beta  = cmplx(beta_re(i),  beta_im(i),  ep)
        transa = transas(i); transb = transbs(i)
        if (transa == 'N') then
            ar = m; ac = k
        else
            ar = k; ac = m
        end if
        if (transb == 'N') then
            br = k; bc = n
        else
            br = n; bc = k
        end if

        call gen_distrib_matrix_z(ar, ac, mb, nb, A_loc, A_glob, seed = 12101 + 43 * i)
        call gen_distrib_matrix_z(br, bc, mb, nb, B_loc, B_glob, seed = 12111 + 43 * i)
        call gen_distrib_matrix_z(m,  n,  mb, nb, C_loc, C0,     seed = 12121 + 43 * i)

        call local_desc(desca, ar, ac, mb, nb)
        call local_desc(descb, br, bc, mb, nb)
        call local_desc(descc, m,  n,  mb, nb)

        call target_pzgemm(transa, transb, m, n, k, alpha, &
                           A_loc, 1, 1, desca, B_loc, 1, 1, descb, &
                           beta,  C_loc, 1, 1, descc)
        call gather_matrix_z(m, n, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(C_ref(m, n))
            C_ref = C0
            call zgemm(transa, transb, m, n, k, alpha, A_glob, ar, B_glob, br, &
                       beta, C_ref, m)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 64.0_ep * 8.0_ep * real(k, ep) * target_eps
            write(label, '(a,a,a,a,a,i0,a,i0,a,i0)') &
                'ta=', transa, ',tb=', transb, &
                ',m=', m, ',n=', n, ',k=', k
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end if
        deallocate(A_loc, B_loc, C_loc, A_glob, B_glob, C0)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgemm
