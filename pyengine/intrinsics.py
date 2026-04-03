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
    # Double-precision specific → generic
    'DCONJG': ('CONJG', False),
    'DIMAG':  ('AIMAG', False),
    'DABS':   ('ABS',   False),
    'DSQRT':  ('SQRT',  False),
    'DEXP':   ('EXP',   False),
    'DLOG':   ('LOG',   False),
    'DSIN':   ('SIN',   False),
    'DCOS':   ('COS',   False),
    'DTAN':   ('TAN',   False),
    'DASIN':  ('ASIN',  False),
    'DACOS':  ('ACOS',  False),
    'DATAN':  ('ATAN',  False),
    'DATAN2': ('ATAN2', False),
    'DSIGN':  ('SIGN',  False),
    'DMAX1':  ('MAX',   False),
    'DMIN1':  ('MIN',   False),
    'DMOD':   ('MOD',   False),
    'DNINT':  ('ANINT', False),
    'IDNINT': ('NINT',  False),
    # Single-precision specific → generic
    'ALOG':   ('LOG',   False),
    'ALOG10': ('LOG10', False),
    'AMAX1':  ('MAX',   False),
    'AMIN1':  ('MIN',   False),
    'AMOD':   ('MOD',   False),
    'CABS':   ('ABS',   False),
    'CSQRT':  ('SQRT',  False),
    'CEXP':   ('EXP',   False),
    'CLOG':   ('LOG',   False),
    # Type conversion — needs KIND arg to preserve target precision
    'DBLE':   ('REAL',  True),
    'DCMPLX': ('CMPLX', True),
    'SNGL':   ('REAL',  True),
    'DFLOAT': ('REAL',  True),
    'FLOAT':  ('REAL',  True),
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
