! The one abort body every Fortran-side libmpiseq stub shares.
!
! The stub files here (mpiseq_qx_stubs.f, mpiseq_ep_stubs.f,
! mpiseq_mw_stubs.f90, mpiseq_ptscotch_stubs.f90) exist only so a
! fully-sequential link resolves; single-rank operation never reaches
! any of them. Each one used to carry its own copy of upstream libseq's
! two-line body — print "should not be called", STOP — which is 44
! copies of one behavior. They call this instead, and the message text
! is written once. The C-side twin is mpiseq_should_not_be_called() in
! mpiseq_c_stubs.c.
!
! Its own translation unit, deliberately: folding it into one of the
! stub files would make every caller pull that file's archive member
! too, and mpiseq_qx_stubs.f's members (PQGETRF, ...) collide with the
! migrated ScaLAPACK archives that legitimately define those names.
! This object defines nothing but the helper.
!
! Free-form (.f90) so both the fixed- and free-form stub files can call
! it; it takes no derived types, so it needs no module and stays a plain
! external subroutine.

subroutine mpiseq_stop_stub(name)
   implicit none
   character(len=*), intent(in) :: name
   write(*,*) 'Error. ' // name // ' should not be called.'
   stop
end subroutine mpiseq_stop_stub
