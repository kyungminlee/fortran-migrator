! ICNTL(18)=3 — user-supplied distributed assembled input for ZMUMPS,
! single rank. Mirror of test_dmumps_dist_input.

program test_zmumps_dist_input
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, nz
    complex(ep), allocatable :: A(:,:), x_true(:), b(:), x_solve(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc)       :: id
    real(ep)                 :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_zmumps_dist_input', target_name)

    call gen_dense_problem_z(n, A, x_true, b, seed = 25001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
    call target_xmumps(id)
    id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    id%ICNTL(18) = 3

    id%N        = n
    id%NZ_loc   = nz
    id%NNZ_loc  = int(nz, kind=8)
    allocate(id%IRN_loc(nz));  id%IRN_loc = irn
    allocate(id%JCN_loc(nz));  id%JCN_loc = jcn
    allocate(id%A_loc(nz));    id%A_loc   = A_trip
    allocate(id%RHS(n));       id%RHS = b

    id%JOB = 6
    call target_xmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(18)=3 failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if

    allocate(x_solve(n));  x_solve = id%RHS
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('icntl18=3', err, tol)

    deallocate(id%IRN_loc, id%JCN_loc, id%A_loc, id%RHS, x_solve)
    nullify(id%IRN_loc, id%JCN_loc, id%A_loc, id%RHS)
    id%JOB = -2;  call target_xmumps(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_zmumps_dist_input
