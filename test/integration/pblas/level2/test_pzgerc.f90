program test_pzgerc
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zgerc
    use pblas_grid,    only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib, only: gen_distrib_matrix_z, gen_distrib_vector_z, &
                             gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pzgerc
    implicit none

    integer, parameter :: ms(*) = [24, 64, 100]
    integer, parameter :: ns(*) = [32, 48, 120]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n
    integer :: desca(9), descx(9), descy(9)
    complex(ep), allocatable :: A_loc(:,:), x_loc(:), y_loc(:)
    complex(ep), allocatable :: A_glob(:,:), x_glob(:), y_glob(:), &
                                A_got(:,:), A_ref(:,:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzgerc', target_name, my_rank)

    alpha = cmplx(0.4_ep, 0.2_ep, ep)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_distrib_matrix_z(m, n, mb, nb, A_loc, A_glob, seed = 5501 + 23 * i)
        call gen_distrib_vector_z(m, mb, x_loc, x_glob, seed = 5511 + 23 * i)
        call gen_distrib_vector_z(n, nb, y_loc, y_glob, seed = 5521 + 23 * i)

        call local_desc(desca, m, n, mb, nb)
        call local_desc(descx, m, 1, mb, 1)
        call local_desc(descy, n, 1, nb, 1)

        call target_pzgerc(m, n, alpha, x_loc, 1, 1, descx, 1, &
                           y_loc, 1, 1, descy, 1, A_loc, 1, 1, desca)
        call gather_matrix_z(m, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(m, n))
            A_ref = A_glob
            call zgerc(m, n, alpha, x_glob, 1, y_glob, 1, A_ref, m)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 32.0_ep * 8.0_ep * target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end if
        deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgerc
