program test_pdsymm
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dsymm
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix, gather_matrix
    use target_pblas,  only: target_name, target_eps, target_pdsymm
    implicit none

    integer, parameter :: ms(*) = [32, 80, 160]
    integer, parameter :: ns(*) = [32, 60, 120]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, info
    integer :: locm_a, locn_a, locm_b, locn_b, locm_c, locn_c
    integer :: lld_a, lld_b, lld_c
    integer :: desca(9), descb(9), descc(9)
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B_glob(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdsymm', target_name, my_rank)

    alpha = 0.7_ep; beta = 0.3_ep
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        ! Side='L': A is m×m symmetric, B is m×n, C is m×n.
        call gen_distrib_matrix(m, m, mb, mb, A_loc, A_glob, seed = 8101 + 29 * i)
        call gen_distrib_matrix(m, n, mb, nb, B_loc, B_glob, seed = 8111 + 29 * i)
        call gen_distrib_matrix(m, n, mb, nb, C_loc, C0,     seed = 8121 + 29 * i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(m, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_b = numroc_local(n, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
        locm_c = locm_b; locn_c = locn_b; lld_c = lld_b

        call descinit_local(desca, m, m, mb, mb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, m, n, mb, nb, 0, 0, my_context, lld_b, info)
        call descinit_local(descc, m, n, mb, nb, 0, 0, my_context, lld_c, info)

        call target_pdsymm('L', 'U', m, n, alpha, &
                           A_loc, 1, 1, desca, B_loc, 1, 1, descb, &
                           beta,  C_loc, 1, 1, descc)
        call gather_matrix(m, n, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(C_ref(m, n))
            C_ref = C0
            call dsymm('L', 'U', m, n, alpha, A_glob, m, B_glob, m, &
                       beta, C_ref, m)
            err = max_rel_err_mat(C_got, C_ref)
            tol = 32.0_ep * 2.0_ep * real(m, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end if
        deallocate(A_loc, B_loc, C_loc, A_glob, B_glob, C0)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdsymm
