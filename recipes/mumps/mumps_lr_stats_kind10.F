C
C  Hand-written extended-precision bridge for MUMPS_LR_STATS — kind10.
C
C  Sibling of mumps_lr_stats_ep.F (kind16). Casts the migrated
C  REAL(KIND=10) input back to DP before calling the upstream
C  MUMPS_LR_STATS implementation.
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
        REAL(KIND=10), INTENT(IN)  :: FLOP_NUMBER
        INTEGER(8), INTENT(OUT) ::
     &                  NB_ENTRIES_FACTOR_withLR
        CALL COMPUTE_GLOBAL_GAINS_DP(NB_ENTRIES_FACTOR,
     &              dble(FLOP_NUMBER), NB_ENTRIES_FACTOR_withLR,
     &              PROKG, MPG)
      END SUBROUTINE COMPUTE_GLOBAL_GAINS
      END MODULE MUMPS_LR_STATS_EP
