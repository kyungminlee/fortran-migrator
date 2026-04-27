! Minimal Fortran shims for the three ScaLAPACK TOOLS integer helpers
! that PBBLAS' pb{q,x}tran / pb{q,x}trnv reference: NUMROC, ICEIL, ILCM.
!
! These live in external/scalapack-2.2.3/TOOLS/{numroc,iceil,ilcm}.f
! upstream, but since the test subtree is forbidden to edit external/
! and only depends on libqpbblas / libqblacs / libqblas (no
! libqscalapack), we re-implement them here byte-for-byte. The
! routines are precision-independent (integer-only), so a single
! definition serves every target. Implemented as plain external
! subprograms (not inside a module) so the linker resolves them via
! the trailing-underscore Fortran ABI used by pbblas.

INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NB, IPROC, ISRCPROC, NPROCS
    INTEGER :: MYDIST, NBLOCKS, EXTRABLKS

    MYDIST  = MOD( NPROCS + IPROC - ISRCPROC, NPROCS )
    NBLOCKS = N / NB
    NUMROC  = ( NBLOCKS / NPROCS ) * NB
    EXTRABLKS = MOD( NBLOCKS, NPROCS )
    IF ( MYDIST .LT. EXTRABLKS ) THEN
        NUMROC = NUMROC + NB
    ELSE IF ( MYDIST .EQ. EXTRABLKS ) THEN
        NUMROC = NUMROC + MOD( N, NB )
    END IF
END FUNCTION NUMROC

INTEGER FUNCTION ICEIL( INUM, IDENOM )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: INUM, IDENOM
    IF ( IDENOM .EQ. 0 ) THEN
        ICEIL = 0
    ELSE
        ICEIL = ( INUM + IDENOM - 1 ) / IDENOM
    END IF
END FUNCTION ICEIL

INTEGER FUNCTION ILCM( M, N )
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M, N
    INTEGER :: A, B, T
    ! Compute GCD(M,N) by Euclid; LCM = M*N/GCD.
    A = ABS(M); B = ABS(N)
    IF ( A .EQ. 0 .OR. B .EQ. 0 ) THEN
        ILCM = 0
        RETURN
    END IF
    DO WHILE ( B .NE. 0 )
        T = MOD( A, B )
        A = B
        B = T
    END DO
    ILCM = ( ABS(M) / A ) * ABS(N)
END FUNCTION ILCM
