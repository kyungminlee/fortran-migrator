program test_pdposv
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dposv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, local_desc
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix, &
                                scatter_matrix
    use target_scalapack, only: target_name, target_eps, target_pdposv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 3
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, info_ref, k
    integer :: desca(9), descb(9)
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B_glob(:,:), B_got(:,:)
    real(ep), allocatable :: A_sym(:,:), A_ref(:,:), B_ref(:,:)
    real(ep) :: err, tol

    call grid_init()
    call report_init('pdposv', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n,    mb, nb, A_loc, A_glob, seed = 12001 + 31*i)
        call gen_distrib_matrix(n, nrhs, mb, nb, B_loc, B_glob, seed = 12011 + 31*i)

        ! Symmetrize and add a diagonal boost so A is SPD.
        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))
        do k = 1, n
            A_sym(k, k) = A_sym(k, k) + real(2 * n, ep)
        end do

        call scatter_matrix(n, n, mb, nb, A_sym, A_loc)
        call local_desc(desca, n, n,    mb, nb)
        call local_desc(descb, n, nrhs, mb, nb)

        call target_pdposv(uplos(i), n, nrhs, A_loc, 1, 1, desca, &
                           B_loc, 1, 1, descb, info)
        call gather_matrix(n, nrhs, mb, nb, B_loc, B_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), B_ref(n, nrhs))
            A_ref = A_sym
            B_ref = B_glob
            call dposv(uplos(i), n, nrhs, A_ref, n, B_ref, n, info_ref)
            err = max_rel_err_mat(B_got, B_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            block
                character(len=48) :: label
                write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
                    ',n=', n, ',nrhs=', nrhs
                call report_case(trim(label), err, tol)
            end block
            deallocate(A_ref, B_ref, B_got)
        end if
        deallocate(A_loc, B_loc, A_glob, B_glob, A_sym)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdposv
