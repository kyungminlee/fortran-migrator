"""Fortran type-specific intrinsic function mappings.

These are universal across all Fortran numerical libraries (BLAS, LAPACK,
ScaLAPACK, MUMPS, etc.) — they are part of the Fortran language, not
library-specific.
"""

# Intrinsic function call replacements.
# Key: old name (type-specific)
# Value: (new_name, needs_kind_arg)
#   needs_kind_arg=False → simple rename: DCONJG(x) → CONJG(x)
#   needs_kind_arg=True  → currently unused; in BLAS/LAPACK context the
#     generic form works because all floating types are at the target KIND.
INTRINSIC_MAP: dict[str, tuple[str, bool]] = {
    # Conjugate, imaginary part, abs, math functions — use generic
    'DCONJG': ('CONJG', False),
    'DIMAG':  ('AIMAG', False),
    'DABS':   ('ABS',   False),
    'DSQRT':  ('SQRT',  False),
    'DEXP':   ('EXP',   False),
    'DLOG':   ('LOG',   False),
    'DSIN':   ('SIN',   False),
    'DCOS':   ('COS',   False),
    'DSIGN':  ('SIGN',  False),
    'DMAX1':  ('MAX',   False),
    'DMIN1':  ('MIN',   False),
    'DNINT':  ('ANINT', False),
    'IDNINT': ('NINT',  False),
    'CABS':   ('ABS',   False),
    # Type conversion — generic form works when all types are at target KIND
    'DBLE':   ('REAL',  False),
    'DCMPLX': ('CMPLX', False),
    'SNGL':   ('REAL',  False),
}

# Replacements for names in INTRINSIC declaration statements.
INTRINSIC_DECL_MAP: dict[str, str] = {
    'DCONJG': 'CONJG',
    'DIMAG':  'AIMAG',
    'DABS':   'ABS',
    'DSQRT':  'SQRT',
    'DBLE':   'REAL',
    'DCMPLX': 'CMPLX',
    'SNGL':   'REAL',
    'CABS':   'ABS',
}
