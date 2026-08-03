C  Per-precision ScaLAPACK stub forwarders for libmpiseq.
C
C  Upstream's extern/MUMPS_5.8.2/libseq/mpi.f provides single-process
C  "should not be called" stubs only for the D-prefixed ScaLAPACK
C  symbols (PDGETRF / PDGETRS / PDPOTRF / PDPOTRS / PDTRTRS). The
C  migrated qmumps / xmumps (kind16) and emumps / ymumps (kind10)
C  archives request the Q/X/E/Y-prefixed equivalents — those resolve
C  through ${LIB_PREFIX}scalapack at the standard mpiexec path, but
C  for a fully-sequential libmpiseq link those symbols have to live
C  somewhere too. This file ships them here in src/ (the
C  extern/ tree is read-only) and is appended to the mpiseq target
C  in src/mpiseq/CMakeLists.txt.
C
C  Despite the _qx_ filename, this file carries all four kind-based
C  prefixes: Q/X (kind16, REAL*16) and E/Y (kind10, REAL*10).
C
C  Each stub mirrors upstream's pattern exactly: print "should not be
C  called", STOP. Single-rank operation never reaches them; they exist
C  only so the link resolves. That body lives once, in
C  mpiseq_stop_stub.f90, and every stub here calls it.
C
C  Multifloats prefixes (M, W) live in mpiseq_mw_stubs.f90 — kept in a
C  separate free-form file because the stubs require TYPE(real64x2) /
C  TYPE(cmplx64x2) from the multifloats Fortran module.
C
C***********************************************************************
      SUBROUTINE PQGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      REAL*16            A( * )
        CALL MPISEQ_STOP_STUB('PQGETRF')
      END SUBROUTINE PQGETRF
C***********************************************************************
      SUBROUTINE PQGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      REAL*16            A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PQGETRS')
      END SUBROUTINE PQGETRS
C***********************************************************************
      SUBROUTINE PQPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      REAL*16            A( * )
        CALL MPISEQ_STOP_STUB('PQPOTRF')
      END SUBROUTINE PQPOTRF
C***********************************************************************
      SUBROUTINE PQPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL*16            A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PQPOTRS')
      END SUBROUTINE PQPOTRS
C***********************************************************************
      SUBROUTINE PQTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL*16            A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PQTRTRS')
      END SUBROUTINE PQTRTRS
C***********************************************************************
      SUBROUTINE PXGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX*32         A( * )
        CALL MPISEQ_STOP_STUB('PXGETRF')
      END SUBROUTINE PXGETRF
C***********************************************************************
      SUBROUTINE PXGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX*32         A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PXGETRS')
      END SUBROUTINE PXGETRS
C***********************************************************************
      SUBROUTINE PXPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      COMPLEX*32         A( * )
        CALL MPISEQ_STOP_STUB('PXPOTRF')
      END SUBROUTINE PXPOTRF
C***********************************************************************
      SUBROUTINE PXPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*32         A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PXPOTRS')
      END SUBROUTINE PXPOTRS
C***********************************************************************
      SUBROUTINE PXTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*32         A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PXTRTRS')
      END SUBROUTINE PXTRTRS
C***********************************************************************
      SUBROUTINE PEGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      REAL*10            A( * )
        CALL MPISEQ_STOP_STUB('PEGETRF')
      END SUBROUTINE PEGETRF
C***********************************************************************
      SUBROUTINE PEGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      REAL*10            A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PEGETRS')
      END SUBROUTINE PEGETRS
C***********************************************************************
      SUBROUTINE PEPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      REAL*10            A( * )
        CALL MPISEQ_STOP_STUB('PEPOTRF')
      END SUBROUTINE PEPOTRF
C***********************************************************************
      SUBROUTINE PEPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL*10            A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PEPOTRS')
      END SUBROUTINE PEPOTRS
C***********************************************************************
      SUBROUTINE PETRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL*10            A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PETRTRS')
      END SUBROUTINE PETRTRS
C***********************************************************************
      SUBROUTINE PYGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX*20         A( * )
        CALL MPISEQ_STOP_STUB('PYGETRF')
      END SUBROUTINE PYGETRF
C***********************************************************************
      SUBROUTINE PYGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX*20         A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PYGETRS')
      END SUBROUTINE PYGETRS
C***********************************************************************
      SUBROUTINE PYPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      COMPLEX*20         A( * )
        CALL MPISEQ_STOP_STUB('PYPOTRF')
      END SUBROUTINE PYPOTRF
C***********************************************************************
      SUBROUTINE PYPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*20         A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PYPOTRS')
      END SUBROUTINE PYPOTRS
C***********************************************************************
      SUBROUTINE PYTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*20         A( * ), B( * )
        CALL MPISEQ_STOP_STUB('PYTRTRS')
      END SUBROUTINE PYTRTRS
