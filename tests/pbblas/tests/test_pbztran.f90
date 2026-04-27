! Test PBZTRAN — distributed conjugate-transpose of a column block to
! a row block. Same scaffolding as test_pbdtran but with TRANS = 'C'
! so the operation is C := A**H + beta*C.
program test_pbztran
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_mat_z
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_matrix_complex
    use pblas_grid,         only: grid_init, grid_exit, my_rank, my_context, &
                                  my_nprow, my_npcol, my_row, my_col, &
                                  numroc_local
    use target_pbblas,      only: target_name, target_eps, target_pbztran
    use mpi
    implicit none

    integer, parameter :: nb = 4
    integer, parameter :: m  = 8
    integer, parameter :: n  = nb
    integer, parameter :: iarow = 0, iacol = 0, icrow = 0, iccol = 0
    integer :: mp, mq, lda, ldc, lwork
    integer :: bytes_per_elem, ig, jg, ierr
    complex(ep), allocatable :: A_glob(:,:), A_loc(:,:), C_loc(:,:), C0_loc(:,:)
    complex(ep), allocatable :: C_glob_got(:,:), C_glob_ref(:,:), C0_glob(:,:)
    complex(ep) :: beta
    real(ep) :: err, tol
    character(len=64) :: label
    integer :: locm_C, locn_C, src_rank, peer
    integer, parameter :: tag_C = 7401
    complex(ep), allocatable :: buf2d(:,:)
    integer :: pc

    call grid_init()
    call report_init('pbztran', target_name, my_rank)

    mp = numroc_local(m, nb, my_row, iarow, my_nprow)
    mq = numroc_local(m, nb, my_col, iccol, my_npcol)
    lda = max(1, mp)
    ldc = max(1, n)
    lwork = max(64, n * m)

    call gen_matrix_complex(m, n, A_glob, seed = 4501)
    allocate(A_loc(lda, max(1, n)))
    A_loc = (0.0_ep, 0.0_ep)
    if (my_col == iacol) then
        do ig = 1, m
            block
                integer :: zero_based, blk, owner_r, il
                zero_based = ig - 1
                blk = zero_based / nb
                owner_r = mod(blk + iarow, my_nprow)
                if (owner_r == my_row) then
                    il = (blk / my_nprow) * nb + mod(zero_based, nb) + 1
                    A_loc(il, 1:n) = A_glob(ig, 1:n)
                end if
            end block
        end do
    end if

    call gen_matrix_complex(n, m, C0_glob, seed = 4601)
    allocate(C_loc(ldc, max(1, mq)), C0_loc(ldc, max(1, mq)))
    C_loc  = (0.0_ep, 0.0_ep)
    C0_loc = (0.0_ep, 0.0_ep)
    if (my_row == icrow .and. mq > 0) then
        do jg = 1, m
            block
                integer :: zero_based, blk, owner_c, jl
                zero_based = jg - 1
                blk = zero_based / nb
                owner_c = mod(blk + iccol, my_npcol)
                if (owner_c == my_col) then
                    jl = (blk / my_npcol) * nb + mod(zero_based, nb) + 1
                    C_loc (1:n, jl) = C0_glob(1:n, jg)
                    C0_loc(1:n, jl) = C0_glob(1:n, jg)
                end if
            end block
        end do
    end if

    beta = (0.5_ep, -0.25_ep)
    call target_pbztran(my_context, 'C', 'C', m, n, nb, A_loc, lda, beta, &
                        C_loc, ldc, iarow, iacol, icrow, iccol, &
                        ncolsA = n, ncolsC = max(1, mq), lwork = lwork)

    bytes_per_elem = storage_size((0.0_ep, 0.0_ep)) / 8
    if (my_rank == 0) then
        allocate(C_glob_got(n, m))
        C_glob_got = (0.0_ep, 0.0_ep)
        do pc = 0, my_npcol - 1
            locm_C = n
            locn_C = numroc_local(m, nb, pc, iccol, my_npcol)
            if (locn_C == 0) cycle
            allocate(buf2d(locm_C, locn_C))
            if (pc == my_col .and. my_row == icrow) then
                buf2d = C_loc(1:locm_C, 1:locn_C)
            else
                src_rank = icrow * my_npcol + pc
                call mpi_recv(buf2d, locm_C * locn_C * bytes_per_elem, &
                              mpi_byte, src_rank, tag_C, mpi_comm_world, &
                              mpi_status_ignore, ierr)
            end if
            do jg = 1, m
                block
                    integer :: zero_based, blk, owner_c, jl
                    zero_based = jg - 1
                    blk = zero_based / nb
                    owner_c = mod(blk + iccol, my_npcol)
                    if (owner_c == pc) then
                        jl = (blk / my_npcol) * nb + mod(zero_based, nb) + 1
                        C_glob_got(1:n, jg) = buf2d(1:n, jl)
                    end if
                end block
            end do
            deallocate(buf2d)
        end do
    else
        if (my_row == icrow .and. mq > 0) then
            peer = 0
            call mpi_send(C_loc, n * mq * bytes_per_elem, mpi_byte, peer, &
                          tag_C, mpi_comm_world, ierr)
        end if
    end if

    if (my_rank == 0) then
        allocate(C_glob_ref(n, m))
        C_glob_ref = transpose(conjg(A_glob)) + beta * C0_glob
        err = max_rel_err_mat_z(C_glob_got, C_glob_ref)
        tol = 64.0_ep * real(m * n, ep) * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(C_glob_got, C_glob_ref)
    end if

    deallocate(A_glob, A_loc, C_loc, C0_loc, C0_glob)

    call report_finalize()
    call grid_exit()
end program test_pbztran
