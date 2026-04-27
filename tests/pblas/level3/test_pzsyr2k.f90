program test_pzsyr2k
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zsyr2k
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pzsyr2k
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: ks(*) = [20, 40, 72]
    integer, parameter :: mb = 8
    integer :: i, n, k, info
    integer :: locm_a, locn_a, locm_b, locn_b, locm_c, locn_c
    integer :: lld_a, lld_b, lld_c
    integer :: desca(9), descb(9), descc(9)
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), B_glob(:,:), C0(:,:), &
                                C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzsyr2k', target_name, my_rank)

    alpha = cmplx(0.5_ep, 0.2_ep, ep); beta = cmplx(0.25_ep, -0.1_ep, ep)
    do i = 1, size(ns)
        n = ns(i); k = ks(i)
        call gen_distrib_matrix_z(n, k, mb, mb, A_loc, A_glob, seed = 14401 + 53 * i)
        call gen_distrib_matrix_z(n, k, mb, mb, B_loc, B_glob, seed = 14411 + 53 * i)
        call gen_distrib_matrix_z(n, n, mb, mb, C_loc, C0,     seed = 14421 + 53 * i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(k, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = locm_a; locn_b = locn_a; lld_b = lld_a
        locm_c = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_c = numroc_local(n, mb, my_col, 0, my_npcol); lld_c = max(1, locm_c)

        call descinit_local(desca, n, k, mb, mb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, n, k, mb, mb, 0, 0, my_context, lld_b, info)
        call descinit_local(descc, n, n, mb, mb, 0, 0, my_context, lld_c, info)

        call target_pzsyr2k('U', 'N', n, k, alpha, &
                            A_loc, 1, 1, desca, B_loc, 1, 1, descb, &
                            beta, C_loc, 1, 1, descc)
        call gather_matrix_z(n, n, mb, mb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(C_ref(n, n))
            C_ref = C0
            call zsyr2k('U', 'N', n, k, alpha, A_glob, n, B_glob, n, &
                        beta, C_ref, n)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 64.0_ep * 8.0_ep * real(k, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',k=', k
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end if
        deallocate(A_loc, B_loc, C_loc, A_glob, B_glob, C0)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzsyr2k
