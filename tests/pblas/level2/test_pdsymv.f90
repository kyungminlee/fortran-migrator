program test_pdsymv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use ref_quad_blas, only: dsymv
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix, gen_distrib_vector, &
                             gather_vector
    use target_pblas,  only: target_name, target_eps, target_pdsymv
    implicit none

    integer, parameter :: ns(*) = [32, 80, 160]
    integer, parameter :: mb = 8
    integer :: i, n, info
    integer :: locm_a, locn_a, locn_x, lld_a, lld_x
    integer :: desca(9), descx(9), descy(9)
    real(ep), allocatable :: A_loc(:,:), x_loc(:), y_loc(:)
    real(ep), allocatable :: A_glob(:,:), x_glob(:), y_glob(:), y_got(:), y_ref(:)
    real(ep) :: alpha, beta, err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pdsymv', target_name, my_rank)

    alpha = 0.6_ep; beta = 0.2_ep
    do i = 1, size(ns)
        n = ns(i)
        ! A is symmetric: generate a general n×n matrix and let PBLAS
        ! read the UPPER half. The reference dsymv also reads UPPER.
        call gen_distrib_matrix(n, n, mb, mb, A_loc, A_glob, seed = 2101 + 17 * i)
        call gen_distrib_vector(n, mb, x_loc, x_glob, seed = 2111 + 17 * i)
        call gen_distrib_vector(n, mb, y_loc, y_glob, seed = 2121 + 17 * i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locn_x = numroc_local(n, mb, my_row, 0, my_nprow); lld_x = max(1, locn_x)

        call descinit_local(desca, n, n, mb, mb, 0, 0, my_context, lld_a, info)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld_x, info)
        call descinit_local(descy, n, 1, mb, 1, 0, 0, my_context, lld_x, info)

        call target_pdsymv('U', n, alpha, A_loc, 1, 1, desca, &
                           x_loc, 1, 1, descx, 1, beta, &
                           y_loc, 1, 1, descy, 1)
        call gather_vector(n, mb, y_loc, y_got)

        if (my_rank == 0) then
            allocate(y_ref(n))
            y_ref = y_glob
            call dsymv('U', n, alpha, A_glob, n, x_glob, 1, beta, y_ref, 1)
            err = max_rel_err_vec(y_got, y_ref)
            tol = 32.0_ep * 2.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(y_ref, y_got)
        end if
        deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdsymv
