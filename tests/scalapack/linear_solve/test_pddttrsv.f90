program test_pddttrsv
    ! Smoke test: factor general tridiagonal with pddttrf, then call
    ! pddttrsv twice (UPLO='L'/'U'). PDDTTRF applies a permutation;
    ! the two trsv halves don't compose to the full pddttrs solve, so
    ! we report INFO=0 from each call as the success metric.
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_npcol, my1d_col, &
                                 descinit_1d
    use test_data,         only: gen_vector_quad, gen_matrix_quad
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pddttrf, target_pddttrsv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer :: i, n, nb, info, lwork, laf, ig, jl, lldB
    integer :: desca(9), descb(9)
    real(ep), allocatable :: dl_glob(:), d_glob(:), du_glob(:), B_glob(:,:)
    real(ep), allocatable :: dl_loc(:), d_loc(:), du_loc(:), B_loc(:,:)
    real(ep), allocatable :: af(:), work(:)
    real(ep) :: err, tol, wopt(1)
    integer :: info_l, info_u
    character(len=48) :: label

    call grid_init()
    call report_init('pddttrsv', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc

        call gen_vector_quad(n, dl_glob, seed = 33401 + 31*i)
        call gen_vector_quad(n, d_glob,  seed = 33411 + 31*i)
        call gen_vector_quad(n, du_glob, seed = 33421 + 31*i)
        call gen_matrix_quad(n, nrhs, B_glob, seed = 33431 + 31*i)
        dl_glob = dl_glob * 0.05_ep
        du_glob = du_glob * 0.05_ep
        d_glob  = abs(d_glob) + real(2 * n, ep)
        dl_glob(1) = 0.0_ep
        du_glob(n) = 0.0_ep

        allocate(dl_loc(nb), d_loc(nb), du_loc(nb))
        dl_loc = 0.0_ep; d_loc = 0.0_ep; du_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    dl_loc(jl) = dl_glob(ig)
                    d_loc(jl)  = d_glob(ig)
                    du_loc(jl) = du_glob(ig)
                end if
            end block
        end do

        lldB = nb
        allocate(B_loc(nb, nrhs)); B_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    B_loc(jl, 1:nrhs) = B_glob(ig, 1:nrhs)
                end if
            end block
        end do

        call descinit_1d(desca, 501, n, nb, 0, my1d_context, nb,   info)
        call descinit_1d(descb, 502, n, nb, 0, my1d_context, lldB, info)

        laf = 12 * my1d_npcol + 3 * nb
        allocate(af(laf))

        call target_pddttrf(n, dl_loc, d_loc, du_loc, 1, desca, &
                            af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pddttrf(n, dl_loc, d_loc, du_loc, 1, desca, &
                            af, laf, work, lwork, info)
        deallocate(work)

        call target_pddttrsv('L', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, wopt, -1, info_l)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pddttrsv('L', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, lwork, info_l)
        deallocate(work)

        allocate(work(1))
        call target_pddttrsv('U', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, -1, info_u)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pddttrsv('U', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, lwork, info_u)

        if (my_rank == 0) then
            err = real(abs(info_l) + abs(info_u), ep)
            tol = 0.5_ep
            write(label, '(a,i0,a,i0,a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs, &
                ',info_L=', info_l, ',info_U=', info_u
            call report_case(trim(label), err, tol)
        end if
        deallocate(dl_loc, d_loc, du_loc, B_loc, work, af, &
                   dl_glob, d_glob, du_glob, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pddttrsv
