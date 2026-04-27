program test_pdtrmv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dtrmv
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix, gen_distrib_vector, &
                             gather_vector
    use target_pblas,  only: target_name, target_eps, target_pdtrmv
    implicit none

    integer, parameter :: ns(*) = [32, 80, 160]
    integer, parameter :: mb = 8
    integer :: i, n, info
    integer :: locm_a, locn_a, locn_x, lld_a, lld_x
    integer :: desca(9), descx(9)
    real(ep), allocatable :: A_loc(:,:), x_loc(:)
    real(ep), allocatable :: A_glob(:,:), x_glob(:), x_got(:), x_ref(:)
    real(ep) :: err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pdtrmv', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, mb, A_loc, A_glob, seed = 5101 + 19 * i)
        call gen_distrib_vector(n, mb, x_loc, x_glob, seed = 5111 + 19 * i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locn_x = numroc_local(n, mb, my_row, 0, my_nprow); lld_x = max(1, locn_x)

        call descinit_local(desca, n, n, mb, mb, 0, 0, my_context, lld_a, info)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld_x, info)

        call target_pdtrmv('U', 'N', 'N', n, A_loc, 1, 1, desca, &
                           x_loc, 1, 1, descx, 1)
        call gather_vector(n, mb, x_loc, x_got)

        if (my_rank == 0) then
            allocate(x_ref(n))
            x_ref = x_glob
            call dtrmv('U', 'N', 'N', n, A_glob, n, x_ref, 1)
            err = max_rel_err_vec(x_got, x_ref)
            tol = 32.0_ep * 2.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(x_ref, x_got)
        end if
        deallocate(A_loc, x_loc, A_glob, x_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdtrmv
