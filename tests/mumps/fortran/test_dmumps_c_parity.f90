! Fortran-vs-C parity test for the QMUMPS path. Same problem driven
! through both target_qmumps (Fortran wrapper, calls QMUMPS which
! calls qmumps_f77_) and the C bridge (qmumps_c which extracts
! struct fields and calls qmumps_f77_). If the bridge correctly
! mirrors the Fortran derived-type extraction, the two solutions
! must be bit-identical.

program test_dmumps_c_parity
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, target_eps, &
                                     dmumps_struc, target_qmumps
    use mpi
    use iso_c_binding,         only: c_int, c_int64_t
    implicit none

    interface
        function c_qmumps_solve(n, nnz, irn, jcn, a_vals, rhs) result(code) &
            bind(C, name='c_qmumps_solve')
            import :: c_int, c_int64_t
            integer(c_int)              :: code
            integer(c_int)              :: n
            integer(c_int64_t)          :: nnz
            integer(c_int)              :: irn(*), jcn(*)
            real(16)                    :: a_vals(*), rhs(*)  ! __float128 layout
        end function c_qmumps_solve
    end interface

    integer, parameter :: n = 16
    integer            :: ierr, nz, c_code
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    real(ep), allocatable :: x_fortran(:), x_c(:)
    type(dmumps_struc)    :: id
    integer(c_int)        :: n_c
    integer(c_int64_t)    :: nnz_c
    real(ep)              :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_dmumps_c_parity', target_name)

    call gen_dense_problem(n, A, x_true, b, seed = 31001)
    call dense_to_triplet (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    ! Path A — Fortran: target_qmumps → QMUMPS → qmumps_f77_
    id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
    call target_qmumps(id)
    id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    id%N   = n
    id%NNZ = int(nz, kind=8)
    allocate(id%IRN(nz));  id%IRN = irn
    allocate(id%JCN(nz));  id%JCN = jcn
    allocate(id%A(nz));    id%A   = A_trip
    allocate(id%RHS(n));   id%RHS = b
    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'Fortran path failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    allocate(x_fortran(n));  x_fortran = id%RHS
    deallocate(id%IRN, id%JCN, id%A, id%RHS)
    nullify(id%IRN, id%JCN, id%A, id%RHS)
    id%JOB = -2;  call target_qmumps(id)

    ! Path B — C bridge: c_qmumps_solve → qmumps_c → qmumps_f77_
    n_c   = n
    nnz_c = int(nz, kind=8)
    allocate(x_c(n));  x_c = b
    c_code = c_qmumps_solve(n_c, nnz_c, irn, jcn, A_trip, x_c)
    if (c_code < 0) then
        write(*, '(a,i0)') 'C path failed, code=', c_code
        error stop 1
    end if

    ! The two paths must be bit-identical. Both feed qmumps_f77_ the
    ! same scalar/array arguments (the bridge just copies fields out
    ! of the struct), so any non-zero diff would mean the bridge
    ! mishandles a field — silent ABI corruption that the per-side
    ! basic tests can't catch.
    err = max_rel_err_vec(x_fortran, x_c)
    call report_case('fortran-vs-c', err, 0.0_ep)

    err = max_rel_err_vec(x_fortran, x_true)
    call report_case('fortran-vs-truth', err, tol)
    err = max_rel_err_vec(x_c, x_true)
    call report_case('c-vs-truth', err, tol)

    deallocate(A, x_true, b, irn, jcn, A_trip, x_fortran, x_c)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_dmumps_c_parity
