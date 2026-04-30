! Basic smoke + tolerance baseline for the migrated DMUMPS via the
! kind16 ${LIB_PREFIX}mumps archive (called as `qmumps` here through
! the target_mumps wrapper).
!
! For each n ∈ {8, 32}:
!   1. Generate a quad-precision dense problem A·x_true = b with a
!      diagonally-dominant unsymmetric A (SYM=0).
!   2. Convert A to MUMPS triplet (IRN/JCN/A_trip).
!   3. Solve via the migrated qmumps with JOB=-1 → JOB=6 → JOB=-2.
!   4. Compare the recovered solution to x_true at REAL(KIND=ep).
!
! The test emits per-case JSON via prec_report. Tolerance is the
! repo-standard n³ scaling (matches tests/lapack/linear_solve/test_dgesv.f90).

program test_dmumps_basic
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, target_eps, &
                                     dmumps_struc, target_qmumps
    use mpi
    implicit none

    ! n=1 covers the degenerate single-node elimination tree (no
    ! merges, no panel updates). n=8 / n=32 give realistic small /
    ! medium baselines.
    integer, parameter :: ns(*) = [1, 8, 32]
    integer            :: ierr, i, n, nz
    real(ep), allocatable :: A(:,:), x_true(:), b(:), x_solve(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    type(dmumps_struc)    :: id
    real(ep)              :: err, tol
    character(len=48)     :: label

    call MPI_INIT(ierr)
    call report_init('test_dmumps_basic', target_name)

    do i = 1, size(ns)
        n = ns(i)
        call gen_dense_problem(n, A, x_true, b, seed = 1001 + 31 * i)
        call dense_to_triplet(A, irn, jcn, A_trip, nz)

        ! Initialize MUMPS instance.
        id%COMM = MPI_COMM_WORLD
        id%PAR  = 1
        id%SYM  = 0
        id%JOB  = -1
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0)') 'JOB=-1 failed, INFOG(1)=', id%INFOG(1)
            error stop 1
        end if

        ! Silence MUMPS diagnostic / info output.
        id%ICNTL(1) = -1
        id%ICNTL(2) = -1
        id%ICNTL(3) = -1
        id%ICNTL(4) = 0

        ! Populate problem data on the host (single-rank => host == only rank).
        id%N    = n
        id%NNZ  = int(nz, kind=8)
        allocate(id%IRN(nz));    id%IRN = irn
        allocate(id%JCN(nz));    id%JCN = jcn
        allocate(id%A(nz));      id%A   = A_trip
        allocate(id%RHS(n));     id%RHS = b

        ! Combined analysis + factorization + solve.
        id%JOB = 6
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0)') 'JOB=6 failed, INFOG(1)=', &
                id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if

        ! On exit, id%RHS holds the solution.
        allocate(x_solve(n))
        x_solve = id%RHS

        err = max_rel_err_vec(x_solve, x_true)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)

        ! Cleanup before next iteration.
        deallocate(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_qmumps(id)

        deallocate(A, x_true, b, x_solve, irn, jcn, A_trip)
    end do

    call report_finalize()
    call MPI_FINALIZE(ierr)
end program test_dmumps_basic
