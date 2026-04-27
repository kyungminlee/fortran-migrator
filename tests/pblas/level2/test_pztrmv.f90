program test_pztrmv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: ztrmv
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gen_distrib_vector_z, &
                             gather_vector_z
    use target_pblas,  only: target_name, target_eps, target_pztrmv
    implicit none

    integer, parameter :: ns(*) = [32, 80, 160]
    integer, parameter :: mb = 8
    integer :: i, n, info
    integer :: locm_a, locn_a, locn_x, lld_a, lld_x
    integer :: desca(9), descx(9)
    complex(ep), allocatable :: A_loc(:,:), x_loc(:)
    complex(ep), allocatable :: A_glob(:,:), x_glob(:), x_got(:), x_ref(:)
    real(ep) :: err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pztrmv', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, mb, A_loc, A_glob, seed = 5901 + 19 * i)
        call gen_distrib_vector_z(n, mb, x_loc, x_glob, seed = 5911 + 19 * i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locn_x = numroc_local(n, mb, my_row, 0, my_nprow); lld_x = max(1, locn_x)

        call descinit_local(desca, n, n, mb, mb, 0, 0, my_context, lld_a, info)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld_x, info)

        call target_pztrmv('U', 'N', 'N', n, A_loc, 1, 1, desca, &
                           x_loc, 1, 1, descx, 1)
        call gather_vector_z(n, mb, x_loc, x_got)

        if (my_rank == 0) then
            allocate(x_ref(n))
            x_ref = x_glob
            call ztrmv('U', 'N', 'N', n, A_glob, n, x_ref, 1)
            err = max_rel_err_vec_z(x_got, x_ref)
            tol = 32.0_ep * 8.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(x_ref, x_got)
        end if
        deallocate(A_loc, x_loc, A_glob, x_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pztrmv
