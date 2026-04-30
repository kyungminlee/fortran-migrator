! ICNTL coverage for ZMUMPS — mirror of test_dmumps_icntl_io.

program test_zmumps_icntl_io
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
    integer, parameter :: scaling_modes(*) = [0, 1, 7, 77]
    integer            :: i, sc
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    complex(ep), allocatable :: x_solve(:)
    type(zmumps_struc)       :: id
    real(ep)                 :: err, tol
    character(len=48)        :: label

    call MPI_INIT(ierr)
    call report_init('test_zmumps_icntl_io', target_name)
    call gen_dense_problem_z(n, A, x_true, b, seed = 26001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    do i = 1, size(scaling_modes)
        sc = scaling_modes(i)
        call init_id(id)
        id%ICNTL(8) = sc
        call attach_dense(id, n, nz, irn, jcn, A_trip, b)
        id%JOB = 6
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0)') 'ICNTL(8)=', sc, &
                ' failed: INFOG(1)=', id%INFOG(1)
            error stop 1
        end if
        allocate(x_solve(n));  x_solve = id%RHS
        err = max_rel_err_vec_z(x_solve, x_true)
        write(label, '(a,i0)') 'icntl8=', sc
        call report_case(trim(label), err, tol)
        call end_id(id);  deallocate(x_solve)
    end do

    call init_id(id)
    id%ICNTL(21) = 1
    call attach_dense(id, n, nz, irn, jcn, A_trip, b)
    id%LSOL_loc = n
    allocate(id%SOL_loc(n))
    allocate(id%ISOL_loc(n))
    id%JOB = 6
    call target_xmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(21)=1 failed: INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    block
        complex(ep) :: x_perm(n)
        integer  :: idx
        do idx = 1, n
            x_perm(id%ISOL_loc(idx)) = id%SOL_loc(idx)
        end do
        err = max_rel_err_vec_z(x_perm, x_true)
    end block
    call report_case('icntl21=1', err, tol)
    deallocate(id%SOL_loc, id%ISOL_loc)
    call end_id(id)

    call init_id(id)
    id%ICNTL(20) = 1
    id%N    = n
    id%NNZ  = int(nz, kind=8)
    allocate(id%IRN(nz));  id%IRN = irn
    allocate(id%JCN(nz));  id%JCN = jcn
    allocate(id%A(nz));    id%A   = A_trip
    id%NRHS    = 1
    id%LRHS    = n
    id%NZ_RHS  = n
    allocate(id%RHS_SPARSE(n));   id%RHS_SPARSE  = b
    allocate(id%IRHS_SPARSE(n))
    allocate(id%IRHS_PTR(2))
    do i = 1, n;  id%IRHS_SPARSE(i) = i;  end do
    id%IRHS_PTR(1) = 1
    id%IRHS_PTR(2) = n + 1
    allocate(id%RHS(n));  id%RHS = cmplx(0.0_ep, 0.0_ep, kind=ep)
    id%JOB = 6
    call target_xmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(20)=1 failed: INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    allocate(x_solve(n));  x_solve = id%RHS
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('icntl20=1', err, tol)
    deallocate(id%RHS_SPARSE, id%IRHS_SPARSE, id%IRHS_PTR, x_solve)
    call end_id(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine init_id(id)
        type(zmumps_struc), intent(inout) :: id
        id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
        call target_xmumps(id)
        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    end subroutine init_id

    subroutine attach_dense(id, n, nz, irn, jcn, A_trip, b)
        type(zmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: n, nz, irn(:), jcn(:)
        complex(ep),        intent(in)    :: A_trip(:), b(:)
        id%N   = n
        id%NNZ = int(nz, kind=8)
        allocate(id%IRN(nz));  id%IRN = irn
        allocate(id%JCN(nz));  id%JCN = jcn
        allocate(id%A(nz));    id%A   = A_trip
        allocate(id%RHS(n));   id%RHS = b
    end subroutine attach_dense

    subroutine end_id(id)
        type(zmumps_struc), intent(inout) :: id
        if (associated(id%IRN)) deallocate(id%IRN)
        if (associated(id%JCN)) deallocate(id%JCN)
        if (associated(id%A))   deallocate(id%A)
        if (associated(id%RHS)) deallocate(id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_xmumps(id)
    end subroutine end_id

end program test_zmumps_icntl_io
