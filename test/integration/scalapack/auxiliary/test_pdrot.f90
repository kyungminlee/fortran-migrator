program test_pdrot
    ! Givens rotation applied to two distributed vectors. Quad-computed
    ! analytic result is the reference (no LAPACK serial proxy needed).
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib,     only: gen_distrib_matrix, gather_matrix
    use target_scalapack,  only: target_name, target_eps, target_pdrot
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, lwork
    integer :: descx(9), descy(9)
    real(ep), allocatable :: x_loc(:,:), x_glob(:,:), x_got(:,:)
    real(ep), allocatable :: y_loc(:,:), y_glob(:,:), y_got(:,:)
    real(ep), allocatable :: x_ref(:), y_ref(:), work(:), wopt(:)
    real(ep) :: cs, sn, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdrot', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        cs = 0.6_ep + 0.01_ep * real(i, ep)
        sn = sqrt(1.0_ep - cs * cs)

        call gen_distrib_matrix(n, 1, mb, nb, x_loc, x_glob, seed = 23601 + 31*i)
        call gen_distrib_matrix(n, 1, mb, nb, y_loc, y_glob, seed = 23611 + 31*i)

        call local_desc(descx, n, 1, mb, nb)
        call local_desc(descy, n, 1, mb, nb)

        allocate(wopt(1))
        call target_pdrot(n, x_loc, 1, 1, descx, 1, y_loc, 1, 1, descy, 1, &
                          cs, sn, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        deallocate(wopt); allocate(work(lwork))
        call target_pdrot(n, x_loc, 1, 1, descx, 1, y_loc, 1, 1, descy, 1, &
                          cs, sn, work, lwork, info)
        deallocate(work)

        call gather_matrix(n, 1, mb, nb, x_loc, x_got)
        call gather_matrix(n, 1, mb, nb, y_loc, y_got)

        if (my_rank == 0) then
            allocate(x_ref(n), y_ref(n))
            x_ref =  cs * x_glob(:, 1) + sn * y_glob(:, 1)
            y_ref = -sn * x_glob(:, 1) + cs * y_glob(:, 1)
            err = max(max_rel_err_vec(x_got(:, 1), x_ref), &
                      max_rel_err_vec(y_got(:, 1), y_ref))
            tol = 32.0_ep * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(x_ref, y_ref, x_got, y_got)
        end if
        deallocate(x_loc, x_glob, y_loc, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdrot
