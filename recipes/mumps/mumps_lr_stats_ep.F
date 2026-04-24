C
C  Hand-written extended-precision bridge for MUMPS_LR_STATS.
C
C  The upstream MUMPS_LR_STATS stores flop / memory counters as DP
C  (arithmetic-agnostic) and is kept verbatim via copy_files. A few
C  of its public subroutines take DP values computed from migrated,
C  working-precision struct fields (e.g. id%RINFOG, id%DKEEP). After
C  migration those fields are REAL(KIND=16); passing them directly
C  to the DP parameters fails with a type mismatch.
C
C  This bridge module re-exports MUMPS_LR_STATS verbatim (via an
C  unqualified USE with rename) and shadows the single routine that
C  needs an EP signature, casting REAL(KIND=16) inputs back to DP
C  before calling the upstream implementation. Migrated callers are
C  redirected from MUMPS_LR_STATS to MUMPS_LR_STATS_EP via
C  ``module_renames`` in the recipe.
C
      MODULE MUMPS_LR_STATS_EP
      USE MUMPS_LR_STATS,
     &    COMPUTE_GLOBAL_GAINS_DP => COMPUTE_GLOBAL_GAINS
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE COMPUTE_GLOBAL_GAINS(NB_ENTRIES_FACTOR,
     &                FLOP_NUMBER, NB_ENTRIES_FACTOR_withLR,
     &                PROKG, MPG)
        INTEGER(8), INTENT(IN) :: NB_ENTRIES_FACTOR
        INTEGER, INTENT(IN)    :: MPG
        LOGICAL, INTENT(IN)    :: PROKG
        REAL(KIND=16), INTENT(IN)  :: FLOP_NUMBER
        INTEGER(8), INTENT(OUT) ::
     &                  NB_ENTRIES_FACTOR_withLR
        CALL COMPUTE_GLOBAL_GAINS_DP(NB_ENTRIES_FACTOR,
     &              dble(FLOP_NUMBER), NB_ENTRIES_FACTOR_withLR,
     &              PROKG, MPG)
      END SUBROUTINE COMPUTE_GLOBAL_GAINS
      END MODULE MUMPS_LR_STATS_EP
