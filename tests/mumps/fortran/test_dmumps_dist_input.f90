! ICNTL(18)=3 — user-supplied distributed assembled matrix entry on
! a single rank.
!
! With ICNTL(18) = 3 MUMPS expects each rank to provide its share of
! the matrix via NZ_loc / NNZ_loc / IRN_loc / JCN_loc / A_loc instead
! of the centralized IRN / JCN / A path. With NPROCS=1 there is only
! one rank, but the dispatcher still routes through the distributed
! code path — exercising the migrator's renaming of the _loc struct
! members and the corresponding F77 formal arguments.
!
! N must still be set on rank 0 (the host), and the analysis still
! happens centrally: ICNTL(18)=3 only changes how the entries are
! supplied during JOB=1, not the ordering / symbolic factorization.

program test_dmumps_dist_input
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, target_eps, &
                                     dmumps_struc, target_qmumps
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, nz
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    real(ep), allocatable :: x_solve(:)
    type(dmumps_struc)    :: id
    real(ep)              :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_dmumps_dist_input', target_name)
    call gen_dense_problem(n, A, x_true, b, seed = 14001)
    call dense_to_triplet (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
    call target_qmumps(id)
    id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    id%ICNTL(18) = 3

    id%N        = n
    id%NZ_loc   = nz
    id%NNZ_loc  = int(nz, kind=8)
    allocate(id%IRN_loc(nz));  id%IRN_loc = irn
    allocate(id%JCN_loc(nz));  id%JCN_loc = jcn
    allocate(id%A_loc(nz));    id%A_loc   = A_trip

    ! RHS still goes through the centralized id%RHS path on the host.
    allocate(id%RHS(n));   id%RHS = b

    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0,a,i0)') 'JOB=6 with ICNTL(18)=3 failed, INFOG(1)=', &
            id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
        error stop 1
    end if

    allocate(x_solve(n));  x_solve = id%RHS
    err = max_rel_err_vec(x_solve, x_true)
    call report_case('icntl18=3', err, tol)

    deallocate(id%IRN_loc, id%JCN_loc, id%A_loc, id%RHS, x_solve)
    id%JOB = -2;  call target_qmumps(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_dmumps_dist_input
