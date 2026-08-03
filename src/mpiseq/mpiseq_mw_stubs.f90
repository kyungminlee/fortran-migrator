! Multifloats (M/W) ScaLAPACK stub forwarders for libmpiseq.
!
! Symmetric to mpiseq_qx_stubs.f for the multifloats prefixes
! M (real, type(real64x2)) and W (complex, type(cmplx64x2)). The
! migrated mmumps / wmumps archives request these P[MW]* symbols at
! link time; under a real-MPI build they resolve through the migrated
! mscalapack archive, but for a fully-sequential libmpiseq link they
! have to live in mpiseq itself. Upstream's extern/MUMPS_5.8.2/
! libseq/mpi.f only ships the PD* / DOUBLE PRECISION forms.
!
! Each stub mirrors upstream's pattern: print "should not be called",
! STOP. Single-rank operation never reaches them; they exist only so
! the link resolves. That body lives once, in mpiseq_stop_stub.f90,
! and every stub here calls it.
!
! Free-form (.f90) so the multifloats module can be USEd for the
! derived-type argument shadows. Built only when NEEDS_MULTIFLOATS is
! on (see the gate in src/mpiseq/CMakeLists.txt).

!-----------------------------------------------------------------------
subroutine PMGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
   use multifloats, only: real64x2
   implicit none
   integer            :: IA, INFO, JA, M, N
   integer            :: DESCA( * ), IPIV( * )
   type(real64x2)     :: A( * )
   call mpiseq_stop_stub('PMGETRF')
end subroutine PMGETRF

!-----------------------------------------------------------------------
subroutine PMGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, &
                    B, IB, JB, DESCB, INFO )
   use multifloats, only: real64x2
   implicit none
   character          :: TRANS
   integer            :: IA, IB, INFO, JA, JB, N, NRHS
   integer            :: DESCA( * ), DESCB( * ), IPIV( * )
   type(real64x2)     :: A( * ), B( * )
   call mpiseq_stop_stub('PMGETRS')
end subroutine PMGETRS

!-----------------------------------------------------------------------
subroutine PMPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
   use multifloats, only: real64x2
   implicit none
   character          :: UPLO
   integer            :: IA, INFO, JA, N
   integer            :: DESCA( * )
   type(real64x2)     :: A( * )
   call mpiseq_stop_stub('PMPOTRF')
end subroutine PMPOTRF

!-----------------------------------------------------------------------
subroutine PMPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA, &
                    B, IB, JB, DESCB, INFO )
   use multifloats, only: real64x2
   implicit none
   character          :: UPLO
   integer            :: IA, IB, INFO, JA, JB, N, NRHS
   integer            :: DESCA( * ), DESCB( * )
   type(real64x2)     :: A( * ), B( * )
   call mpiseq_stop_stub('PMPOTRS')
end subroutine PMPOTRS

!-----------------------------------------------------------------------
subroutine PMTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA, &
                    B, IB, JB, DESCB, INFO )
   use multifloats, only: real64x2
   implicit none
   character          :: DIAG, TRANS, UPLO
   integer            :: IA, IB, INFO, JA, JB, N, NRHS
   integer            :: DESCA( * ), DESCB( * )
   type(real64x2)     :: A( * ), B( * )
   call mpiseq_stop_stub('PMTRTRS')
end subroutine PMTRTRS

!-----------------------------------------------------------------------
subroutine PWGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
   use multifloats, only: cmplx64x2
   implicit none
   integer            :: IA, INFO, JA, M, N
   integer            :: DESCA( * ), IPIV( * )
   type(cmplx64x2)    :: A( * )
   call mpiseq_stop_stub('PWGETRF')
end subroutine PWGETRF

!-----------------------------------------------------------------------
subroutine PWGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, &
                    B, IB, JB, DESCB, INFO )
   use multifloats, only: cmplx64x2
   implicit none
   character          :: TRANS
   integer            :: IA, IB, INFO, JA, JB, N, NRHS
   integer            :: DESCA( * ), DESCB( * ), IPIV( * )
   type(cmplx64x2)    :: A( * ), B( * )
   call mpiseq_stop_stub('PWGETRS')
end subroutine PWGETRS

!-----------------------------------------------------------------------
subroutine PWPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
   use multifloats, only: cmplx64x2
   implicit none
   character          :: UPLO
   integer            :: IA, INFO, JA, N
   integer            :: DESCA( * )
   type(cmplx64x2)    :: A( * )
   call mpiseq_stop_stub('PWPOTRF')
end subroutine PWPOTRF

!-----------------------------------------------------------------------
subroutine PWPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA, &
                    B, IB, JB, DESCB, INFO )
   use multifloats, only: cmplx64x2
   implicit none
   character          :: UPLO
   integer            :: IA, IB, INFO, JA, JB, N, NRHS
   integer            :: DESCA( * ), DESCB( * )
   type(cmplx64x2)    :: A( * ), B( * )
   call mpiseq_stop_stub('PWPOTRS')
end subroutine PWPOTRS

!-----------------------------------------------------------------------
subroutine PWTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA, &
                    B, IB, JB, DESCB, INFO )
   use multifloats, only: cmplx64x2
   implicit none
   character          :: DIAG, TRANS, UPLO
   integer            :: IA, IB, INFO, JA, JB, N, NRHS
   integer            :: DESCA( * ), DESCB( * )
   type(cmplx64x2)    :: A( * ), B( * )
   call mpiseq_stop_stub('PWTRTRS')
end subroutine PWTRTRS
