! Test BLACS_PINFO and BLACS_GET direct probes.
!
! BLACS_PINFO(MYPNUM, NPROCS) returns this rank's MPI rank and the
! total nproc; subsequent calls must return identical values
! (idempotent).
!
! BLACS_GET(ICONTXT=-1, WHAT=0, VAL) retrieves the default system
! context derived from MPI_COMM_WORLD; calling it with WHAT=10
! retrieves the topology for broadcasts (a small integer).
!
! This program runs after grid_init() has already configured a grid;
! the probes here exercise the bare BLACS-API contract independently
! of the precision-prefixed reduction / broadcast tests.
program test_blacs_pinfo
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my_context
    use target_blacs,      only: target_name
    use mpi
    implicit none

    interface
        subroutine blacs_pinfo(mypnum, nprocs)
            integer, intent(out) :: mypnum, nprocs
        end subroutine
        subroutine blacs_get(icontxt, what, val)
            integer, intent(in)  :: icontxt, what
            integer, intent(out) :: val
        end subroutine
    end interface

    integer :: rk1, rk2, np1, np2, ctxt_default, topo
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)

    call grid_init()
    call report_init('blacs_pinfo', target_name, my_rank)

    ! ── Idempotency: two consecutive PINFO calls must agree, and
    !    must agree with what grid_init recorded.
    call blacs_pinfo(rk1, np1)
    call blacs_pinfo(rk2, np2)
    err = 0.0_ep
    if (rk1 /= rk2 .or. np1 /= np2) err = BAD
    if (rk1 /= my_rank .or. np1 /= my_nproc) err = BAD
    if (my_rank == 0) call report_case('pinfo_idempotent', err, tol)

    ! ── BLACS_GET WHAT=0: default system context. Should return a
    !    valid handle that is *different* from -1 (-1 is the
    !    sentinel meaning "use default") and consistent across calls.
    call blacs_get(-1, 0, ctxt_default)
    err = 0.0_ep
    if (ctxt_default == -1) err = BAD
    if (my_rank == 0) call report_case('get_default_context', err, tol)

    ! ── BLACS_GET WHAT=10: broadcast topology. Returns the default
    !    topology character code (as INTEGER). Must be representable.
    call blacs_get(my_context, 10, topo)
    err = 0.0_ep
    ! topo is a printable ASCII character code; accept any 0..127.
    if (topo < 0 .or. topo > 127) err = BAD
    if (my_rank == 0) call report_case('get_broadcast_topology', err, tol)

    call report_finalize()
    call grid_exit()
end program test_blacs_pinfo
