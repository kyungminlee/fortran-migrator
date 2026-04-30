C
C  Hand-written extended-precision reallocator module — kind10 target.
C
C  Sibling of mumps_memory_mod_ep.F (kind16). Bridges callers migrated
C  to extended-precision arithmetic (E / Y for kind10) to the
C  reallocator contract provided by the upstream MUMPS_MEMORY_MOD.
C
C  See mumps_memory_mod_ep.F for the rationale; only the EP variants
C  differ — kind10 uses REAL(KIND=10) / COMPLEX(KIND=10).
C
      MODULE MUMPS_MEMORY_MOD_EP
      USE MUMPS_MEMORY_MOD, ONLY:
     &     MUMPS_IREALLOC,  MUMPS_IREALLOC8,
     &     MUMPS_I8REALLOC, MUMPS_I8REALLOC8,
     &     MUMPS_SREALLOC,  MUMPS_DREALLOC,
     &     MUMPS_CREALLOC,  MUMPS_ZREALLOC,
     &     MUMPS_IDEALLOC,  MUMPS_I8DEALLOC,
     &     MUMPS_MEMORY_SET_DATA_SIZES
      IMPLICIT NONE
      INTERFACE MUMPS_REALLOC
      MODULE PROCEDURE MUMPS_IREALLOC
      MODULE PROCEDURE MUMPS_SREALLOC, MUMPS_DREALLOC
      MODULE PROCEDURE MUMPS_CREALLOC, MUMPS_ZREALLOC
      MODULE PROCEDURE MUMPS_EREALLOC, MUMPS_YREALLOC
      END INTERFACE
      INTERFACE MUMPS_DEALLOC
      MODULE PROCEDURE MUMPS_IDEALLOC
      END INTERFACE
      CONTAINS
      SUBROUTINE MUMPS_EREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      REAL(KIND=10), POINTER       :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      REAL(KIND=10), POINTER       :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      INTEGER(8)                   :: ELSIZE
      CHARACTER(len=60)            :: FMTA, FMTD
      ELSIZE = int(storage_size(1.0_10),8)/8_8
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*ELSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ELSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ELSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+int(MINSIZE,8)*ELSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_EREALLOC
      SUBROUTINE MUMPS_YREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      COMPLEX(KIND=10), POINTER    :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      COMPLEX(KIND=10), POINTER    :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      INTEGER(8)                   :: ELSIZE
      CHARACTER(len=60)            :: FMTA, FMTD
      ELSIZE = int(storage_size((1.0_10,0.0_10)),8)/8_8
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*ELSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ELSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ELSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+int(MINSIZE,8)*ELSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_YREALLOC
      END MODULE MUMPS_MEMORY_MOD_EP
