program test_pzdttrsv
    ! Smoke test: complex general tridiag pzdttrf + pzdttrsv halves.
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_npcol, my1d_col, &
                                 descinit_1d
    use test_data,         only: gen_vector_quad, gen_matrix_quad
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pzdttrf, target_pzdttrsv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer :: i, n, nb, info, lwork, laf, ig, jl, lldB
    integer :: desca(9), descb(9)
    real(ep),    allocatable :: re1(:), im1(:), re2(:), im2(:), re3(:), im3(:)
    real(ep),    allocatable :: B_re(:,:), B_im(:,:)
    complex(ep), allocatable :: dl_glob(:), d_glob(:), du_glob(:), B_glob(:,:)
    complex(ep), allocatable :: dl_loc(:), d_loc(:), du_loc(:), B_loc(:,:)
    complex(ep), allocatable :: af(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    integer :: info_l, info_u
    character(len=48) :: label

    call grid_init()
    call report_init('pzdttrsv', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc

        call gen_vector_quad(n, re1, seed = 33801 + 31*i)
        call gen_vector_quad(n, im1, seed = 33811 + 31*i)
        call gen_vector_quad(n, re2, seed = 33821 + 31*i)
        call gen_vector_quad(n, im2, seed = 33831 + 31*i)
        call gen_vector_quad(n, re3, seed = 33841 + 31*i)
        call gen_vector_quad(n, im3, seed = 33851 + 31*i)
        call gen_matrix_quad(n, nrhs, B_re, seed = 33861 + 31*i)
        call gen_matrix_quad(n, nrhs, B_im, seed = 33871 + 31*i)
        allocate(dl_glob(n), d_glob(n), du_glob(n), B_glob(n, nrhs))
        dl_glob = cmplx(re1, im1, ep) * 0.05_ep
        du_glob = cmplx(re3, im3, ep) * 0.05_ep
        d_glob  = cmplx(abs(re2) + real(2 * n, ep), 0.0_ep, ep)
        dl_glob(1) = (0.0_ep, 0.0_ep)
        du_glob(n) = (0.0_ep, 0.0_ep)
        B_glob = cmplx(B_re, B_im, ep)

        allocate(dl_loc(nb), d_loc(nb), du_loc(nb))
        dl_loc = (0.0_ep, 0.0_ep); d_loc = (0.0_ep, 0.0_ep); du_loc = (0.0_ep, 0.0_ep)
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
        allocate(B_loc(nb, nrhs)); B_loc = (0.0_ep, 0.0_ep)
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

        call target_pzdttrf(n, dl_loc, d_loc, du_loc, 1, desca, &
                            af, laf, wopt, -1, info)
        lwork = max(1, int(real(wopt(1))))
        allocate(work(lwork))
        call target_pzdttrf(n, dl_loc, d_loc, du_loc, 1, desca, &
                            af, laf, work, lwork, info)
        deallocate(work)

        call target_pzdttrsv('L', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, wopt, -1, info_l)
        lwork = max(1, int(real(wopt(1))))
        allocate(work(lwork))
        call target_pzdttrsv('L', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, lwork, info_l)
        deallocate(work)

        allocate(work(1))
        call target_pzdttrsv('U', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, -1, info_u)
        lwork = max(1, int(real(work(1))))
        deallocate(work); allocate(work(lwork))
        call target_pzdttrsv('U', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, lwork, info_u)

        if (my_rank == 0) then
            err = real(abs(info_l) + abs(info_u), ep)
            tol = 0.5_ep
            write(label, '(a,i0,a,i0,a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs, &
                ',info_L=', info_l, ',info_U=', info_u
            call report_case(trim(label), err, tol)
        end if
        deallocate(dl_loc, d_loc, du_loc, B_loc, work, af, &
                   re1, im1, re2, im2, re3, im3, &
                   dl_glob, d_glob, du_glob, &
                   B_re, B_im, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzdttrsv
