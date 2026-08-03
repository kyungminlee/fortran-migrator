program test_pdtran
    ! pdtran: C := beta*C + alpha*A^T. Reference is element-wise transpose.
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,    only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib, only: gen_distrib_matrix, gather_matrix
    use target_pblas,  only: target_name, target_eps, target_pdtran
    implicit none

    integer, parameter :: ms(*) = [32, 64, 96]
    integer, parameter :: ns(*) = [40, 48, 80]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n
    integer :: desca(9), descc(9)
    real(ep), allocatable :: A_loc(:,:), C_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdtran', target_name, my_rank)

    alpha = 0.7_ep; beta = 0.3_ep
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        ! C is m×n; A is n×m (since result uses A^T).
        call gen_distrib_matrix(n, m, mb, nb, A_loc, A_glob, seed = 17101 + 29*i)
        call gen_distrib_matrix(m, n, mb, nb, C_loc, C0,    seed = 17111 + 29*i)

        call local_desc(desca, n, m, mb, nb)
        call local_desc(descc, m, n, mb, nb)

        call target_pdtran(m, n, alpha, A_loc, 1, 1, desca, &
                           beta, C_loc, 1, 1, descc)
        call gather_matrix(m, n, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(C_ref(m, n))
            C_ref = beta * C0 + alpha * transpose(A_glob)
            err = max_rel_err_mat(C_got, C_ref)
            tol = 16.0_ep * real(max(m, n), ep) * target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end if
        deallocate(A_loc, C_loc, A_glob, C0)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdtran
