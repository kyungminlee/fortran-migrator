program test_pdpttrf
    ! PD tridiagonal factor (pdpttrf) — compare modified d, e against
    ! LAPACK dpttrf. Sibling of test_pdpttrs which exercises factor +
    ! solve; this driver isolates the factor so RESULT.md gets a
    ! pdpttrf-keyed report.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dpttrf
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_col, descinit_1d
    use test_data,         only: gen_vector_quad
    use mpi
    use target_scalapack,  only: target_name, target_eps, target_pdpttrf
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer :: i, n, nb, info, info_ref, lwork, laf, ig, jl
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer, parameter :: tag = 7551
    integer :: desca(9)
    real(ep), allocatable :: d_glob(:), e_glob(:)
    real(ep), allocatable :: d_loc(:),  e_loc(:)
    real(ep), allocatable :: af(:), work(:)
    real(ep), allocatable :: d_got(:), e_got(:), buf_d(:), buf_e(:)
    real(ep), allocatable :: d_ref(:), e_ref(:)
    real(ep) :: err_d, err_e, err, tol, wopt(1)
    character(len=48) :: label

    call grid_init()
    call report_init('pdpttrf', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc

        call gen_vector_quad(n, d_glob, seed = 19501 + 31*i)
        call gen_vector_quad(n, e_glob, seed = 19511 + 31*i)
        d_glob       = abs(d_glob) + real(2 * n, ep)   ! diag-dominant PD
        e_glob       = e_glob * 0.05_ep
        e_glob(n)    = 0.0_ep

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

        call descinit_1d(desca, 501, n, nb, 0, my1d_context, nb, info)

        laf = 12 * my_nproc + 3 * nb
        allocate(af(laf))

        call target_pdpttrf(n, d_loc, e_loc, 1, desca, &
                            af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpttrf(n, d_loc, e_loc, 1, desca, &
                            af, laf, work, lwork, info)
        deallocate(work)

        bytes_per_elem = storage_size(0.0_ep) / 8
        if (my_rank == 0) then
            allocate(d_got(n), e_got(n))
            d_got = 0.0_ep; e_got = 0.0_ep
            do src_rank = 0, my_nproc - 1
                allocate(buf_d(nb), buf_e(nb))
                if (src_rank == my_rank) then
                    buf_d = d_loc; buf_e = e_loc
                else
                    call mpi_recv(buf_d, nb * bytes_per_elem, mpi_byte, &
                                  src_rank, tag,     mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                    call mpi_recv(buf_e, nb * bytes_per_elem, mpi_byte, &
                                  src_rank, tag + 1, mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                end if
                do ig = 1, n
                    block
                        integer :: blk, owner
                        blk = (ig - 1) / nb; owner = blk
                        if (owner == src_rank) then
                            jl = mod(ig - 1, nb) + 1
                            d_got(ig) = buf_d(jl)
                            e_got(ig) = buf_e(jl)
                        end if
                    end block
                end do
                deallocate(buf_d, buf_e)
            end do
        else
            peer = 0
            call mpi_send(d_loc, nb * bytes_per_elem, mpi_byte, &
                          peer, tag,     mpi_comm_world, ierr)
            call mpi_send(e_loc, nb * bytes_per_elem, mpi_byte, &
                          peer, tag + 1, mpi_comm_world, ierr)
        end if

        if (my_rank == 0) then
            allocate(d_ref(n), e_ref(n))
            d_ref = d_glob
            e_ref(1:n-1) = e_glob(1:n-1)
            e_ref(n)     = 0.0_ep
            call dpttrf(n, d_ref, e_ref, info_ref)
            err_d = max_rel_err_vec(d_got(1:n),     d_ref(1:n))
            err_e = max_rel_err_vec(e_got(1:n-1),   e_ref(1:n-1))
            err   = max(err_d, err_e)
            tol   = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(d_ref, e_ref, d_got, e_got)
        end if
        deallocate(d_loc, e_loc, af, d_glob, e_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdpttrf
