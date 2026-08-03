program test_pzdotu
    use prec_kinds,    only: ep
    use compare,       only: rel_err_scalar_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zdotu
    use pblas_grid,    only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib, only: gen_distrib_vector_z
    use target_pblas,  only: target_name, target_eps, target_pzdotu
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n
    integer :: descx(9), descy(9)
    complex(ep), allocatable :: x_loc(:), y_loc(:), x_glob(:), y_glob(:)
    complex(ep) :: ref, got
    real(ep) :: err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pzdotu', target_name, my_rank)

    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector_z(n, mb, x_loc, x_glob, seed = 2801 + 7 * i)
        call gen_distrib_vector_z(n, mb, y_loc, y_glob, seed = 2901 + 7 * i)

        call local_desc(descx, n, 1, mb, 1)
        call local_desc(descy, n, 1, mb, 1)

        got = (0.0_ep, 0.0_ep)
        call target_pzdotu(n, got, x_loc, 1, 1, descx, 1, &
                           y_loc, 1, 1, descy, 1)

        if (my_rank == 0) then
            ref = zdotu(n, x_glob, 1, y_glob, 1)
            err = rel_err_scalar_z(got, ref)
            tol = 32.0_ep * 8.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
        end if
        deallocate(x_loc, y_loc, x_glob, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzdotu
