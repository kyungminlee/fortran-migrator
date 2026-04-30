! Fortran-vs-C parity test for the XMUMPS path. Mirror of
! test_dmumps_c_parity.

program test_zmumps_c_parity
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps
    use mpi
    use iso_c_binding,         only: c_int, c_int64_t
    implicit none

    interface
        function c_xmumps_solve(n, nnz, irn, jcn, a_vals, rhs) result(code) &
            bind(C, name='c_xmumps_solve')
            import :: c_int, c_int64_t
            integer(c_int)              :: code
            integer(c_int)              :: n
            integer(c_int64_t)          :: nnz
            integer(c_int)              :: irn(*), jcn(*)
            complex(16)                 :: a_vals(*), rhs(*)
        end function c_xmumps_solve
    end interface

    integer, parameter :: n = 16
    integer            :: ierr, nz, c_code
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    complex(ep), allocatable :: x_fortran(:), x_c(:)
    type(zmumps_struc)       :: id
    integer(c_int)           :: n_c
    integer(c_int64_t)       :: nnz_c
    real(ep)                 :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_zmumps_c_parity', target_name)

    call gen_dense_problem_z(n, A, x_true, b, seed = 32001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
    call target_xmumps(id)
    id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    id%N   = n
    id%NNZ = int(nz, kind=8)
    allocate(id%IRN(nz));  id%IRN = irn
    allocate(id%JCN(nz));  id%JCN = jcn
    allocate(id%A(nz));    id%A   = A_trip
    allocate(id%RHS(n));   id%RHS = b
    id%JOB = 6
    call target_xmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'Fortran path failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    allocate(x_fortran(n));  x_fortran = id%RHS
    deallocate(id%IRN, id%JCN, id%A, id%RHS)
    nullify(id%IRN, id%JCN, id%A, id%RHS)
    id%JOB = -2;  call target_xmumps(id)

    n_c   = n
    nnz_c = int(nz, kind=8)
    allocate(x_c(n));  x_c = b
    c_code = c_xmumps_solve(n_c, nnz_c, irn, jcn, A_trip, x_c)
    if (c_code < 0) then
        write(*, '(a,i0)') 'C path failed, code=', c_code
        error stop 1
    end if

    err = max_rel_err_vec_z(x_fortran, x_c)
    call report_case('fortran-vs-c', err, 0.0_ep)
    err = max_rel_err_vec_z(x_fortran, x_true)
    call report_case('fortran-vs-truth', err, tol)
    err = max_rel_err_vec_z(x_c, x_true)
    call report_case('c-vs-truth', err, tol)

    deallocate(A, x_true, b, irn, jcn, A_trip, x_fortran, x_c)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_zmumps_c_parity
