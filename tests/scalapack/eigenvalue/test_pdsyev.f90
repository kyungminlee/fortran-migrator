program test_pdsyev
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dsyev
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdsyev
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9), descz(9)
    real(ep), allocatable :: A_loc(:,:), Z_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), A_sym(:,:), A_ref(:,:)
    real(ep), allocatable :: w_got(:), w_ref(:), work(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label
    integer :: ig, jg, owner_r, owner_c, il, jl

    call grid_init()
    call report_init('pdsyev', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 9701 + 31*i)

        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        if (locm_a > 0 .and. locn_a > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_sym(ig, jg)
                end do
            end do
        end if
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descz, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(Z_loc(max(1, locm_a), max(1, locn_a)))
        Z_loc = 0.0_ep
        lwork = max(1, 8 * n + 2 * max(locm_a, locn_a))
        allocate(w_got(n), work(lwork))
        call target_pdsyev('N', 'U', n, A_loc, 1, 1, desca, w_got, &
                           Z_loc, 1, 1, descz, work, lwork, info)

        if (my_rank == 0) then
            allocate(A_ref(n, n), w_ref(n), work_ref(max(1, 64 * n)))
            A_ref = A_sym
            call dsyev('N', 'U', n, A_ref, n, w_ref, work_ref, size(work_ref), info_ref)
            err = max_rel_err_vec(w_got, w_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, w_ref, work_ref)
        end if
        deallocate(A_loc, Z_loc, A_glob, A_sym, w_got, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdsyev
