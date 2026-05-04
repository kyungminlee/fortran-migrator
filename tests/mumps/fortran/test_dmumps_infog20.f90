! Cross-target sanity check on INFOG(20) — peak number of
! real-precision words used during factorization.
!
! Why this exists: each migrator target ships hand-written byte-
! accounting overrides for `mumps_memory_mod` / `mumps_lr_stats`
! (recipe `overrides:`). The override constants encode bytes-per-
! real, bytes-per-complex etc., and getting those wrong is silent —
! the factor still works numerically, but the reported memory /
! work counters drift. INFOG(20) is the most direct user-visible
! signal: it reports the structural element count (in target-
! precision words) and should be roughly target-independent for a
! fixed (n, SYM, JOB) problem because the multifrontal tree shape
! depends on the matrix structure, not the precision.
!
! The bound below is a structural one (n² … 50·n²) rather than a
! kind16 baseline ±5%, so the test can ship before a per-target
! baseline has been captured. Tighten to ±5% against a captured
! reference once enough builds (kind10 + kind16) have logged the
! value into the JSON report.
!
! Multifloats coverage is deferred — tests/mumps/CMakeLists.txt
! gates multifloats out at the executable level
! (_MUMPS_TESTS_SKIP_EXECUTABLES), and Group A in
! tests/mumps/TODO.md tracks the harness rework needed before this
! test can run there. Re-evaluate the bound as a strict-±5%
! cross-target check once Group A lands.

program test_dmumps_infog20
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, &
                                     report_finalize, report_check_status
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, dmumps_struc, target_qmumps, &
                                     q2t_r
    use mpi
    implicit none

    integer, parameter :: n = 32
    integer            :: ierr, nz
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    type(dmumps_struc) :: id
    integer            :: real_words
    real(ep)           :: err, tol
    character(len=64)  :: label

    call MPI_INIT(ierr)
    call report_init('test_dmumps_infog20', target_name)

    call gen_dense_problem(n, A, x_true, b, seed = 4099)
    call dense_to_triplet(A, irn, jcn, A_trip, nz)

    id%COMM = MPI_COMM_WORLD
    id%PAR  = 1
    id%SYM  = 0
    id%JOB  = -1
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'JOB=-1 failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    id%ICNTL(1:3) = -1
    id%ICNTL(4)   = 0

    id%N    = n
    id%NNZ  = int(nz, kind=8)
    allocate(id%IRN(nz));  id%IRN = irn
    allocate(id%JCN(nz));  id%JCN = jcn
    allocate(id%A(nz));    id%A   = q2t_r(A_trip)
    allocate(id%RHS(n));   id%RHS = q2t_r(b)

    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0,a,i0)') 'JOB=6 failed, INFOG(1)=', &
            id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
        error stop 1
    end if

    real_words = id%INFOG(20)

    ! Structural bound: a dense LU on n=32 needs at least ~n² real
    ! words (full L+U) and shouldn't exceed ~50·n² even with MUMPS's
    ! internal padding. Anything outside that range strongly suggests
    ! the target's byte-accounting overrides are mis-sized.
    write(label, '(a,i0,a,i0)') 'infog20=', real_words, ' n=', n
    if (real_words >= n*n .and. real_words <= 50 * n * n) then
        err = 0.0_ep
    else
        err = 1.0_ep
    end if
    tol = 0.5_ep
    call report_case(trim(label), err, tol)

    deallocate(id%IRN, id%JCN, id%A, id%RHS)
    nullify(id%IRN, id%JCN, id%A, id%RHS)
    id%JOB = -2
    call target_qmumps(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)

    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

end program test_dmumps_infog20
