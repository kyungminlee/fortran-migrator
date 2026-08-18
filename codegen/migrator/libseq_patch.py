"""MUMPS libseq ``mpi.f`` patch — a recipe-specific staging hack.

Extracted verbatim from ``__main__.py`` (Cluster 3) as part of the migrator
file-restructuring refactor. Behaviour is unchanged; ``cmd_stage`` imports
``patch_libseq_mpi_f`` from here.
"""
from pathlib import Path

# Derived-type sentinels for the libmpiseq stand-in: ``mpiseq_c_stubs.c``
# encodes a derived MPI datatype handle as this tag OR'd with the type's
# total size in bytes. Keep in sync with ``src/mpiseq/mpiseq_dtype_tag.h``,
# which defines the same encoding for the C side (MPI_Type_c2f /
# MPI_Op_c2f in Intel mpi.h are passthrough casts, so the Fortran handle
# is the same value).
_MPISEQ_DTYPE_TAG = 0x10000000
_FLOAT64X2_BYTES = 16
_COMPLEX64X2_BYTES = 32

# ``MUMPS_COPY_*`` helpers, as (subroutine name, element declaration,
# elements-per-item). The multifloats pair copies through REAL(KIND=8)
# because the staged libseq has no multifloats type: a float64x2 is two
# doubles, a complex64x2 is four. Note the names do NOT track the kinds
# — ``MUMPS_COPY_COMPLEX32`` is COMPLEX(KIND=16) and
# ``MUMPS_COPY_COMPLEX20`` is COMPLEX(KIND=10); the number is the byte
# count, not the kind.
_COPY_HELPERS: tuple[tuple[str, str, int], ...] = (
    ('MUMPS_COPY_REAL16', 'REAL(KIND=16)', 1),
    ('MUMPS_COPY_COMPLEX32', 'COMPLEX(KIND=16)', 1),
    ('MUMPS_COPY_REAL10', 'REAL(KIND=10)', 1),
    ('MUMPS_COPY_COMPLEX20', 'COMPLEX(KIND=10)', 1),
    ('MUMPS_COPY_FLOAT64X2', 'REAL(KIND=8)', 2),
    ('MUMPS_COPY_COMPLEX64X2', 'REAL(KIND=8)', 4),
)


def _copy_helper_text(name: str, decl: str, mult: int) -> str:
    """Emit one ``MUMPS_COPY_*`` subroutine (fixed-form Fortran)."""
    m = '' if mult == 1 else f'{mult}*'
    return (
        f'      SUBROUTINE {name}( S, R, N, SS, RS )\n'
        '      IMPLICIT NONE\n'
        '      INTEGER N, SS, RS\n'
        f'      {decl} S({m}N),R({m}N)\n'
        '      INTEGER I\n'
        f'      DO I = 1, {m}N\n'
        f'        R(I+{m}RS) = S(I+{m}SS)\n'
        '      END DO\n'
        '      RETURN\n'
        f'      END SUBROUTINE {name}\n'
    )


def _patch_mumps_copy_dispatch(src: str) -> str:
    """Extend ``MUMPS_COPY``'s datatype dispatch.

    Upstream's ``MUMPS_COPY`` only knows the standard MPI datatypes; the
    migrated qxmumps archive passes MPI_REAL16 (Intel MPI = 1275072555)
    for kind16 reductions. The new cases are *inserted* immediately
    before the existing ``ELSE / IERR=1`` fallthrough — they must
    precede it, so this is an anchored replace, not an append.
    """
    fallthrough = '      ELSE\n        IERR=1\n        RETURN\n      END IF'
    extra_dispatch = (
        # kind16: REAL(16) / COMPLEX(16)
        '      ELSE IF ( DATATYPE .EQ. MPI_REAL16 ) THEN\n'
        '      CALL MUMPS_COPY_REAL16( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        '      ELSE IF ( DATATYPE .EQ. MPI_COMPLEX32 ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX32( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        # kind10: 80-bit extended real / complex map to MPI's long
        # double tokens (no MPI_REAL10 in standard MPI).
        '      ELSE IF ( DATATYPE .EQ. MPI_LONG_DOUBLE ) THEN\n'
        '      CALL MUMPS_COPY_REAL10( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        '      ELSE IF ( DATATYPE .EQ. MPI_C_LONG_DOUBLE_COMPLEX ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX20( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        # multifloats: derived-type sentinels (see _MPISEQ_DTYPE_TAG).
        f'      ELSE IF ( DATATYPE .EQ. '
        f'{_MPISEQ_DTYPE_TAG | _FLOAT64X2_BYTES} ) THEN\n'
        '      CALL MUMPS_COPY_FLOAT64X2( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        f'      ELSE IF ( DATATYPE .EQ. '
        f'{_MPISEQ_DTYPE_TAG | _COMPLEX64X2_BYTES} ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX64X2( SENDBUF, RECVBUF, CNT, SS, RS )\n'
    )
    if 'MPI_REAL16' not in src and fallthrough in src:
        src = src.replace(fallthrough, extra_dispatch + fallthrough, 1)
    return src


def _append_copy_helpers(src: str) -> str:
    """Append the ``MUMPS_COPY_*`` helpers the new dispatch cases call."""
    if 'SUBROUTINE MUMPS_COPY_REAL16' in src:
        return src
    extra_helpers = '\n' + ''.join(
        _copy_helper_text(*entry) for entry in _COPY_HELPERS)
    return src.rstrip() + '\n' + extra_helpers


def _append_pchk_stubs(src: str) -> str:
    """Append ScaLAPACK descriptor-check stubs.

    Base libseq mpi.f stubs pchk2mat (MUMPS calls it) but NOT pchk1mat /
    globchk. Our migrated *typed* ScaLAPACK archives (ey/qx/mw) are full
    netlib ports that reference all three, co-located in the reference
    ScaLAPACK's single pchkxmat.f.o object. In a seq link, pulling that
    object in for pchk1mat_/globchk_ also drags its pchk2mat_, colliding
    with libseq's pchk2mat_. Stubbing pchk1mat / globchk here makes libseq
    fully cover pchkxmat.f.o, so the reference object is never extracted and
    there is no duplicate symbol (no --allow-multiple-definition). Like
    pchk2mat these are descriptor checks reached only from inside pXgetrf/
    pXpotrf, which MUMPS never calls at np=1 (the root is factored by
    sequential LAPACK) — so the STOP body never executes.
    """
    extra_pchk = """
      SUBROUTINE pchk1mat( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     &                     DESCAPOS0, NEXTRA, EX, EXPOS, INFO )
      IMPLICIT NONE
      INTEGER            DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA,
     &                   NAPOS0, NEXTRA
      INTEGER            DESCA( * ), EX( NEXTRA ), EXPOS( NEXTRA )
        WRITE(*,*) 'Error. PCHK1MAT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pchk1mat
      SUBROUTINE globchk( ICTXT, N, X, LDX, IWORK, INFO )
      IMPLICIT NONE
      INTEGER            ICTXT, INFO, LDX, N
      INTEGER            IWORK( * ), X( LDX, * )
        WRITE(*,*) 'Error. GLOBCHK should not be called.'
        STOP
      RETURN
      END SUBROUTINE globchk
"""
    if 'SUBROUTINE pchk1mat' in src:
        return src
    return src.rstrip() + '\n' + extra_pchk


def patch_libseq_mpi_f(path: Path) -> None:
    """Extend libseq's ``MUMPS_COPY`` with MPI_REAL16 / MPI_COMPLEX32
    cases so reductions on REAL(KIND=16) / COMPLEX(KIND=16) buffers
    dispatch correctly under our libmpiseq variant, and stub the two
    ScaLAPACK descriptor checks libseq is missing.

    Patches the staged copy at ``_mpiseq_src/mpi.f``; upstream's
    ``extern/MUMPS_5.9.1/libseq/mpi.f`` stays read-only. BLACS /
    ScaLAPACK forwarders inside the same file are deliberately KEPT
    — libmpiseq stands in for those archives in the ``_seq`` test
    link, and the real BLACS / ScaLAPACK archives aren't linked there
    so there's no duplicate-symbol collision.
    """
    src = path.read_text()
    src = _patch_mumps_copy_dispatch(src)
    src = _append_copy_helpers(src)
    src = _append_pchk_stubs(src)
    path.write_text(src)
