program test_pzpbtrsv
    ! Smoke test: complex Hermitian-PD banded pzpbtrf + pzpbtrsv halves.
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_col, descinit_1d
    use mpi
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pzpbtrf, target_pzpbtrsv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: bw = 2, lda = bw + 1
    integer :: i, n, nb, info, lwork, laf, ig, jg, jl, lldB
    integer :: desca(9), descb(9)
    complex(ep), allocatable :: A_full(:,:), M(:,:), B_glob(:,:)
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    complex(ep), allocatable :: af(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    integer :: info_u, info_t
    character(len=48) :: label

    call grid_init()
    call report_init('pzpbtrsv', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc
        allocate(A_full(n, n), M(n, n), B_glob(n, nrhs))
        block
            integer :: k, j, s, t
            M = (0.0_ep, 0.0_ep)
            do j = 1, n
                do k = max(1, j - bw), min(n, j + bw)
                    s = 33901 + 31 * i + 7 * j + 11 * k
                    t = 33911 + 31 * i + 7 * j + 11 * k
                    M(k, j) = cmplx(real(mod(s, 1009), ep) / 1009.0_ep, &
                                     real(mod(t, 1019), ep) / 1019.0_ep, ep)
                end do
                do k = 1, nrhs
                    s = 33921 + 31 * i + 7 * j + 11 * k
                    t = 33931 + 31 * i + 7 * j + 11 * k
                    B_glob(j, k) = cmplx(real(mod(s, 1029), ep) / 1029.0_ep, &
                                          real(mod(t, 1031), ep) / 1031.0_ep, ep)
                end do
            end do
        end block
        A_full = matmul(conjg(transpose(M)), M)
        block
            integer :: k, j
            do k = 1, n
                A_full(k, k) = A_full(k, k) + cmplx(real(n, ep), 0.0_ep, ep)
            end do
            do j = 1, n
                do k = 1, n
                    if (abs(k - j) > bw) A_full(k, j) = (0.0_ep, 0.0_ep)
                end do
            end do
        end block

        allocate(A_loc(lda, nb)); A_loc = (0.0_ep, 0.0_ep)
        do jg = 1, n
            block
                integer :: blk, owner, ig2
                blk = (jg - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(jg - 1, nb) + 1
                    do ig2 = max(1, jg - bw), jg
                        A_loc(bw + 1 + ig2 - jg, jl) = A_full(ig2, jg)
                    end do
                end if
            end block
        end do

        lldB = nb
        allocate(B_loc(lldB, nrhs)); B_loc = (0.0_ep, 0.0_ep)
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

        call descinit_1d(desca, 501, n, nb, 0, my1d_context, lda,  info)
        call descinit_1d(descb, 502, n, nb, 0, my1d_context, lldB, info)

        laf = (nb + 2 * bw) * bw + 16
        allocate(af(laf))

        call target_pzpbtrf('U', n, bw, A_loc, 1, desca, &
                            af, laf, wopt, -1, info)
        lwork = max(1, int(real(wopt(1))))
        allocate(work(lwork))
        call target_pzpbtrf('U', n, bw, A_loc, 1, desca, &
                            af, laf, work, lwork, info)
        deallocate(work)

        call target_pzpbtrsv('U', 'N', n, bw, nrhs, A_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, wopt, -1, info_u)
        lwork = max(1, int(real(wopt(1))))
        allocate(work(lwork))
        call target_pzpbtrsv('U', 'N', n, bw, nrhs, A_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, lwork, info_u)
        deallocate(work)

        allocate(work(1))
        call target_pzpbtrsv('U', 'C', n, bw, nrhs, A_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, -1, info_t)
        lwork = max(1, int(real(work(1))))
        deallocate(work); allocate(work(lwork))
        call target_pzpbtrsv('U', 'C', n, bw, nrhs, A_loc, 1, desca, &
                             B_loc, 1, descb, af, laf, work, lwork, info_t)

        if (my_rank == 0) then
            err = real(abs(info_u) + abs(info_t), ep)
            tol = 0.5_ep
            write(label, '(a,i0,a,i0,a,i0,a,i0,a,i0)') 'n=', n, ',bw=', bw, &
                ',nrhs=', nrhs, ',info_N=', info_u, ',info_C=', info_t
            call report_case(trim(label), err, tol)
        end if
        deallocate(A_loc, B_loc, work, af, A_full, M, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzpbtrsv
