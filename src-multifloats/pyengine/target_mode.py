from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class TargetMode:
    """Specifies how the migrator should transform Fortran source code."""

    # Identity
    name: str

    # Type declaration output
    real_type: str
    complex_type: str

    # Literal replacement strategy
    literal_mode: str                  # 'kind_suffix' or 'constructor'
    kind_suffix: Optional[int]         # 10 or 16 for KIND, None for multifloats
    real_constructor: Optional[str]    # None for KIND, 'float64x2' for MF
    complex_constructor: Optional[str] # None for KIND, 'complex128x2' for MF

    # Intrinsic call transformation
    intrinsic_mode: str                # 'add_kind' or 'wrap_constructor'

    # Prefix map for renaming
    prefix_map: dict[str, str]

    # Free-form kind parameter
    kind_param_value: Optional[str]    # '10' or '16' for KIND, None for MF

    # Module / constants support
    module_name: Optional[str]         # None for KIND, 'multifloats' for MF
    known_constants: dict[str, str]    # {'ZERO': 'MF_ZERO', ...}
    la_constants_map: dict[str, str]   # la_constants rename map

    @property
    def is_kind_based(self) -> bool:
        return self.kind_suffix is not None


def kind_target(kind: int) -> TargetMode:
    """Construct TargetMode for standard KIND=10 or KIND=16 compilation."""
    from .prefix_classifier import PREFIX_MAP
    
    if kind not in PREFIX_MAP:
        raise ValueError(f"Unsupported kind: {kind}")
        
    return TargetMode(
        name=f"kind{kind}",
        real_type=f"REAL(KIND={kind})",
        complex_type=f"COMPLEX(KIND={kind})",
        literal_mode="kind_suffix",
        kind_suffix=kind,
        real_constructor=None,
        complex_constructor=None,
        intrinsic_mode="add_kind",
        prefix_map=PREFIX_MAP[kind],
        kind_param_value=str(kind),
        module_name=None,
        known_constants={},
        la_constants_map={}
    )


def multifloats_target(
    module: str = 'multifloats',
    real_type: str = 'float64x2',
    complex_type: str = 'complex128x2',
    prefix_style: str = 'wide',
) -> TargetMode:
    """Construct TargetMode for multifloats library migration."""
    # Prefix convention: W for Real (Wide), U for Complex (Ultra/Wide-Complex)
    if prefix_style == 'wide':
        # Base map from KIND=16, but override R and C prefixes
        # Note: PREFIX_MAP maps logical type to character
        # Standard: 'R' -> 'D', 'C' -> 'Z' for double
        # Quad: 'R' -> 'Q', 'C' -> 'X'
        # Multifloats: 'R' -> 'W', 'C' -> 'U'
        prefix_map = {
            'R': 'W',
            'C': 'U',
            # Int/Logical are unchanged, but we typically don't rename them.
            # Keep defaults just in case
            'I': 'I',
            'L': 'L',
            'S': 'S', # Character
        }
    else:
        raise ValueError(f"Unknown prefix_style: {prefix_style}")

    # Local FP variables in BLAS/LAPACK that are conventionally named
    # after these constants and DATA-initialized to literal values are
    # supplied as named constants by the multifloats module. The
    # migrator filters them from declarations and substitutes references.
    # Names like RTMIN/RTMAX are intentionally NOT included here: in
    # several LAPACK files (e.g. dlartg.f90, dnrm2.f90) they are LOCAL
    # variables that store ``sqrt(safmin)`` / ``sqrt(safmax/2)`` and
    # must remain locally declared and assigned.
    known_constants = {
        'ZERO': 'MF_ZERO',
        'ONE': 'MF_ONE',
        'TWO': 'MF_TWO',
        'HALF': 'MF_HALF',
        'EIGHT': 'MF_EIGHT',
    }

    # Names that LAPACK files import via ``USE LA_CONSTANTS, ONLY: ...``
    # and that have a multifloats equivalent. ``rtmin``/``rtmax`` are
    # NOT renamed here because some LAPACK routines (dlartg, dnrm2)
    # also declare local variables with the same names; mass-renaming
    # them would create assignments to module PARAMETERs.
    la_constants_map = {
        'zero': 'MF_ZERO',
        'half': 'MF_HALF',
        'one': 'MF_ONE',
        'two': 'MF_TWO',
        'eight': 'MF_EIGHT',
        'safmin': 'MF_SAFMIN',
        'safmax': 'MF_SAFMAX',
        'tsml': 'MF_TSML',
        'tbig': 'MF_TBIG',
        'ssml': 'MF_SSML',
        'sbig': 'MF_SBIG',
    }

    return TargetMode(
        name="multifloats",
        real_type=f"TYPE({real_type})",
        complex_type=f"TYPE({complex_type})",
        literal_mode="constructor",
        kind_suffix=None,
        real_constructor=real_type,
        complex_constructor=complex_type,
        intrinsic_mode="wrap_constructor",
        prefix_map=prefix_map,
        kind_param_value=None,
        module_name=module,
        known_constants=known_constants,
        la_constants_map=la_constants_map
    )
