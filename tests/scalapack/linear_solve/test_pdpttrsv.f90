program test_pdpttrsv
    ! Smoke test: factor a SPD tridiagonal with pdpttrf, then call
    ! pdpttrsv twice (UPLO='L' / UPLO='U'). Differential precision
    ! against LAPACK is non-trivial because PDPTTRF applies a
    ! permutation P (P*T*P^T = L*D*L^T) and pdpttrsv applies only
    ! one triangular half — composing the two halves does not
    ! reproduce the full pdpttrs solve. We verify INFO=0 from each
    ! call instead and report it as a binary pass/fail.
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_npcol, my1d_col, &
                                 descinit_1d
    use test_data,         only: gen_vector_quad, gen_matrix_quad
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pdpttrf, target_pdpttrsv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer :: i, n, nb, info, lwork, laf, ig, jl, lldB
    integer :: desca(9), descb(9)
    real(ep), allocatable :: d_glob(:), e_glob(:), B_glob(:,:)
    real(ep), allocatable :: d_loc(:),  e_loc(:),  B_loc(:,:)
    real(ep), allocatable :: af(:), work(:)
    real(ep) :: err, tol, wopt(1)
    integer :: info_l, info_u
    character(len=48) :: label

    call grid_init()
    call report_init('pdpttrsv', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc

        call gen_vector_quad(n, d_glob, seed = 33301 + 31*i)
        call gen_vector_quad(n, e_glob, seed = 33311 + 31*i)
        call gen_matrix_quad(n, nrhs, B_glob, seed = 33321 + 31*i)
        d_glob = abs(d_glob) + real(2 * n, ep)
        e_glob = e_glob * 0.05_ep
        e_glob(n) = 0.0_ep

        allocate(d_loc(nb), e_loc(nb))
        d_loc = 0.0_ep; e_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    d_loc(jl) = d_glob(ig)
                    e_loc(jl) = e_glob(ig)
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

        call target_pdpttrf(n, d_loc, e_loc, 1, desca, &
                            af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpttrf(n, d_loc, e_loc, 1, desca, &
                            af, laf, work, lwork, info)
        deallocate(work)

        call target_pdpttrsv('L', n, nrhs, d_loc, e_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, wopt, -1, info_l)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpttrsv('L', n, nrhs, d_loc, e_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, lwork, info_l)
        deallocate(work)

        allocate(work(1))
        call target_pdpttrsv('U', n, nrhs, d_loc, e_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, -1, info_u)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdpttrsv('U', n, nrhs, d_loc, e_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, lwork, info_u)

        if (my_rank == 0) then
            err = real(abs(info_l) + abs(info_u), ep)
            tol = 0.5_ep
            write(label, '(a,i0,a,i0,a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs, &
                ',info_L=', info_l, ',info_U=', info_u
            call report_case(trim(label), err, tol)
        end if
        deallocate(d_loc, e_loc, B_loc, work, af, &
                   d_glob, e_glob, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdpttrsv
