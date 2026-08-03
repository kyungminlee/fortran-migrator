C     ep_-privatized twins of libseq's BLACS / ScaLAPACK abort stubs.
C
C     Upstream libseq's mpi.f defines pristine-named abort stubs
C     (blacs_gridinit, DESCINIT, ...) for the grid entry points the
C     parallel MUMPS archives reference but a single-rank run never
C     executes (the type-3 parallel root is never selected at
C     NPROCS=1).  The migrated extended-precision MUMPS archives
C     (m/w, q/x, e/y) reach those routines through the privatized
C     ep_ BLACS engine names instead (doc/guide/mkl-coexistence.md),
C     so a libmpiseq link needs the same abort stubs under the ep_
C     spelling.  Exactly the four names extended MUMPS binds
C     (codegen/recipes/mumps.yaml extra_renames).  Always folded into
C     libmpiseq: baseline kind4/kind8 consumers never reference
C     ep_* and the archive member stays unpulled.
C***********************************************************************
      SUBROUTINE ep_blacs_gridinit( CNTXT, C, NPROW, NPCOL )
      IMPLICIT NONE
      INTEGER CNTXT, NPROW, NPCOL
      CHARACTER C
        CALL MPISEQ_STOP_STUB('EP_BLACS_GRIDINIT')
      END SUBROUTINE ep_blacs_gridinit
C***********************************************************************
      SUBROUTINE ep_blacs_gridinfo( CNTXT, NPROW, NPCOL, MYROW, MYCOL )
      IMPLICIT NONE
      INTEGER CNTXT, NPROW, NPCOL, MYROW, MYCOL
        CALL MPISEQ_STOP_STUB('EP_BLACS_GRIDINFO')
      END SUBROUTINE ep_blacs_gridinfo
C***********************************************************************
      SUBROUTINE ep_blacs_gridexit( CNTXT )
      IMPLICIT NONE
      INTEGER CNTXT
        CALL MPISEQ_STOP_STUB('EP_BLACS_GRIDEXIT')
      END SUBROUTINE ep_blacs_gridexit
C***********************************************************************
      SUBROUTINE ep_descinit( DESC, M, N, MB, NB, IRSRC, ICSRC,
     &           ICTXT, LLD, INFO )
      IMPLICIT NONE
      INTEGER ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER DESC( * )
        CALL MPISEQ_STOP_STUB('EP_DESCINIT')
      END SUBROUTINE ep_descinit
