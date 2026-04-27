! Test PBDVECADD — local vector add Y := alpha*op(X) + beta*Y.
! Same pattern as test_pbdmatadd: routine performs no communication,
! each rank uses distinct local data, comparison is local + a global
! pass/fail reduction so all ranks agree on the exit code.
program test_pbdvecadd
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_vec
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_vector_quad
    use pblas_grid,         only: grid_init, grid_exit, my_rank, my_context
    use target_pbblas,      only: target_name, target_eps, target_pbdvecadd
    use mpi
    implicit none

    integer, parameter :: cases(*) = [16, 100, 1000]
    character(len=2), parameter :: modes(*) = ['V ', 'C ']
    integer :: n, i, k, ierr, fail_local, fail_global
    real(ep), allocatable :: X(:), Y0(:), Y_got(:), Y_ref(:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call grid_init()
    call report_init('pbdvecadd', target_name, my_rank)

    alpha = 1.7_ep; beta = -0.5_ep

    do k = 1, size(modes)
        do i = 1, size(cases)
            n = cases(i)
            call gen_vector_quad(n, X,    seed = 711 + 11 * i + 100 * my_rank + 1000 * k)
            call gen_vector_quad(n, Y0,   seed = 811 + 11 * i + 100 * my_rank + 1000 * k)
            allocate(Y_got(n), Y_ref(n))
            Y_got = Y0; Y_ref = Y0

            call target_pbdvecadd(my_context, modes(k)(1:1), n, alpha, X, 1, &
                                  beta, Y_got, 1)

            ! Reference: real-data MODE='C' is identical to vector add.
            Y_ref = alpha * X + beta * Y_ref

            err = max_rel_err_vec(Y_got, Y_ref)
            tol = 64.0_ep * real(n, ep) * target_eps

            fail_local = 0
            if (err > tol) fail_local = 1
            call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                               mpi_max, mpi_comm_world, ierr)

            if (my_rank == 0) then
                write(label, '(a,a,a,i0,a,i0)') 'mode=', trim(modes(k)), &
                    ',n=', n, ',anyfail=', fail_global
                if (fail_global /= 0 .and. err <= tol) err = 2.0_ep * tol
                call report_case(trim(label), err, tol)
            end if
            deallocate(X, Y0, Y_got, Y_ref)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pbdvecadd
