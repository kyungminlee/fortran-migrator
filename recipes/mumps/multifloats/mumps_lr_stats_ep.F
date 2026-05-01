C
C  Hand-written extended-precision bridge for MUMPS_LR_STATS —
C  multifloats target. FLOP_NUMBER is TYPE(real64x2) and is cast to
C  DP via dble() (the multifloats module exports this as a generic
C  to-double conversion) before delegating to the upstream routine.
C
      MODULE MUMPS_LR_STATS_EP
      USE multifloats
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
        TYPE(real64x2), INTENT(IN)  :: FLOP_NUMBER
        INTEGER(8), INTENT(OUT) ::
     &                  NB_ENTRIES_FACTOR_withLR
        CALL COMPUTE_GLOBAL_GAINS_DP(NB_ENTRIES_FACTOR,
     &              dble(FLOP_NUMBER), NB_ENTRIES_FACTOR_withLR,
     &              PROKG, MPG)
      END SUBROUTINE COMPUTE_GLOBAL_GAINS
      END MODULE MUMPS_LR_STATS_EP
