C
C  Hand-written extended-precision bridge for MUMPS_LR_STATS — multifloats.
C
C  Sibling of mumps_lr_stats_ep.F (kind16). Casts the migrated
C  TYPE(real64x2) input back to DP before calling the upstream
C  MUMPS_LR_STATS implementation. The multifloats module supplies a
C  generic ``dble()`` overload that returns the high limb of a
C  double-double — sufficient for accounting / reporting.
C
      MODULE MUMPS_LR_STATS_EP
      USE multifloats, ONLY: real64x2, dble
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
