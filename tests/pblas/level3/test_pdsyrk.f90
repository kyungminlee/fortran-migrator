program test_pdsyrk
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use ref_quad_blas, only: dsyrk
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix, gather_matrix
    use target_pblas,  only: target_name, target_eps, target_pdsyrk
    implicit none

    integer, parameter :: ns(*) = [32, 80, 128]
    integer, parameter :: ks(*) = [24, 48, 100]
    integer, parameter :: mb = 8
    integer :: i, n, k, info
    integer :: locm_a, locn_a, locm_c, locn_c, lld_a, lld_c
    integer :: desca(9), descc(9)
    real(ep), allocatable :: A_loc(:,:), C_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdsyrk', target_name, my_rank)

    alpha = 0.5_ep; beta = 0.25_ep
    do i = 1, size(ns)
        n = ns(i); k = ks(i)
        call gen_distrib_matrix(n, k, mb, mb, A_loc, A_glob, seed = 9201 + 31 * i)
        call gen_distrib_matrix(n, n, mb, mb, C_loc, C0,     seed = 9211 + 31 * i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(k, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_c = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_c = numroc_local(n, mb, my_col, 0, my_npcol); lld_c = max(1, locm_c)

        call descinit_local(desca, n, k, mb, mb, 0, 0, my_context, lld_a, info)
        call descinit_local(descc, n, n, mb, mb, 0, 0, my_context, lld_c, info)

        call target_pdsyrk('U', 'N', n, k, alpha, &
                           A_loc, 1, 1, desca, beta, C_loc, 1, 1, descc)
        call gather_matrix(n, n, mb, mb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(C_ref(n, n))
            C_ref = C0
            call dsyrk('U', 'N', n, k, alpha, A_glob, n, beta, C_ref, n)
            err = max_rel_err_mat(C_got, C_ref)
            tol = 32.0_ep * 2.0_ep * real(k, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',k=', k
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end if
        deallocate(A_loc, C_loc, A_glob, C0)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdsyrk
