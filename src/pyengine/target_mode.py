from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


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
    kind_suffix: Optional[int]         # 10 or 16 for KIND, None for module-based
    real_constructor: Optional[str]    # None for KIND, 'float64x2' for MF
    complex_constructor: Optional[str] # None for KIND, 'complex128x2' for MF

    # Intrinsic call transformation
    intrinsic_mode: str                # 'add_kind' or 'wrap_constructor'

    # Prefix map for renaming
    prefix_map: dict[str, str]

    # Module / constants support
    module_name: Optional[str]         # None for KIND, 'multifloats' for MF
    known_constants: dict[str, str]    # {'ZERO': 'MF_ZERO', ...}
    la_constants_map: dict[str, str]   # la_constants rename map
    la_constants_suffix: str           # '_MF' or '_EP'

    # Module public names (for building USE...ONLY clauses)
    module_type_names: frozenset[str] = field(default_factory=frozenset)
    module_constant_names: frozenset[str] = field(default_factory=frozenset)
    module_generic_names: frozenset[str] = field(default_factory=frozenset)
    module_operator_generics: tuple[str, ...] = ()

    # C-interop fields (None for KIND targets, populated for module-based).
    # Used by c_migrator to substitute types / MPI handles / reduction ops
    # in cloned BLACS and PBLAS C sources.
    c_real_type: Optional[str] = None         # 'float64x2_t'
    c_complex_type: Optional[str] = None      # 'complex128x2_t'
    c_mpi_real: Optional[str] = None          # 'MPI_FLOAT64X2'
    c_mpi_complex: Optional[str] = None       # 'MPI_COMPLEX128X2'
    c_mpi_sum_real: Optional[str] = None      # 'MPI_DD_SUM'
    c_mpi_sum_complex: Optional[str] = None   # 'MPI_ZZ_SUM'
    c_header: Optional[str] = None            # 'multifloats_bridge.h'

    @property
    def is_kind_based(self) -> bool:
        return self.kind_suffix is not None

    @property
    def module_public_names(self) -> frozenset[str]:
        """All public names exported by the module (types + constants + generics)."""
        return self.module_type_names | self.module_constant_names | self.module_generic_names


def load_target(name_or_path: str, project_root: Path | None = None) -> TargetMode:
    """Load a TargetMode from a YAML file or built-in name.

    Resolution order:
      1. If name_or_path is a path to an existing .yaml file, load it directly.
      2. Look for targets/{name_or_path}.yaml relative to project_root.
      3. If name_or_path is a bare integer (e.g. "16"), try targets/kind{n}.yaml.
    """
    if project_root is None:
        # Default: repo root is two levels up from this file
        # (src/pyengine/target_mode.py → repo root)
        project_root = Path(__file__).resolve().parent.parent.parent

    # 1. Direct path
    p = Path(name_or_path)
    if p.suffix in ('.yaml', '.yml') and p.exists():
        return _load_target_yaml(p)

    # 2. Named lookup
    targets_dir = project_root / 'targets'
    candidate = targets_dir / f'{name_or_path}.yaml'
    if candidate.exists():
        return _load_target_yaml(candidate)

    # 3. Bare integer → kindN
    if name_or_path.isdigit():
        candidate = targets_dir / f'kind{name_or_path}.yaml'
        if candidate.exists():
            return _load_target_yaml(candidate)

    raise FileNotFoundError(
        f"No target YAML found for '{name_or_path}'. "
        f"Searched: {targets_dir / f'{name_or_path}.yaml'}"
    )


def _load_target_yaml(path: Path) -> TargetMode:
    """Parse a target YAML file and construct a TargetMode."""
    with open(path) as f:
        d = yaml.safe_load(f)

    # Module section: can be null (KIND targets) or a dict (module-based)
    mod = d.get('module')
    if isinstance(mod, dict):
        module_name = mod['name']
        module_type_names = frozenset(mod.get('type_names') or [])
        module_constant_names = frozenset(mod.get('constant_names') or [])
        module_generic_names = frozenset(mod.get('generic_names') or [])
        module_operator_generics = tuple(mod.get('operator_generics') or [])
    else:
        module_name = mod  # None or a bare string (backward compat)
        module_type_names = frozenset()
        module_constant_names = frozenset()
        module_generic_names = frozenset()
        module_operator_generics = ()

    c = d.get('c_interop') or {}

    return TargetMode(
        name=d['name'],
        real_type=d['types']['real'],
        complex_type=d['types']['complex'],
        literal_mode=d['literals']['mode'],
        kind_suffix=d['literals'].get('kind_suffix'),
        real_constructor=d['literals'].get('real_constructor'),
        complex_constructor=d['literals'].get('complex_constructor'),
        intrinsic_mode=d['intrinsics']['mode'],
        prefix_map=d['prefixes'],
        module_name=module_name,
        known_constants=d.get('known_constants') or {},
        la_constants_map=d.get('la_constants_map') or {},
        la_constants_suffix=d.get('la_constants_suffix', ''),
        module_type_names=module_type_names,
        module_constant_names=module_constant_names,
        module_generic_names=module_generic_names,
        module_operator_generics=module_operator_generics,
        c_real_type=c.get('real_type'),
        c_complex_type=c.get('complex_type'),
        c_mpi_real=c.get('mpi_real'),
        c_mpi_complex=c.get('mpi_complex'),
        c_mpi_sum_real=c.get('mpi_sum_real'),
        c_mpi_sum_complex=c.get('mpi_sum_complex'),
        c_header=c.get('header'),
    )


# -- Deprecated factory functions (kept for backward compatibility) ----------

def kind_target(kind: int) -> TargetMode:
    """Construct TargetMode for standard KIND=10 or KIND=16 compilation.

    .. deprecated:: Use ``load_target('kind10')`` or ``load_target('kind16')`` instead.
    """
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
        module_name=None,
        known_constants={},
        la_constants_map={},
        la_constants_suffix='_EP',
    )


def multifloats_target(
    module: str = 'multifloats',
    real_type: str = 'float64x2',
    complex_type: str = 'complex128x2',
    prefix_style: str = 'wide',
) -> TargetMode:
    """Construct TargetMode for multifloats library migration.

    .. deprecated:: Use ``load_target('multifloats')`` instead.
    """
    if prefix_style == 'wide':
        prefix_map = {
            'R': 'DD',
            'C': 'ZZ',
            'I': 'I',
            'L': 'L',
            'S': 'S',
        }
    else:
        raise ValueError(f"Unknown prefix_style: {prefix_style}")

    known_constants = {
        'ZERO': 'MF_ZERO',
        'ONE': 'MF_ONE',
        'TWO': 'MF_TWO',
        'HALF': 'MF_HALF',
        'EIGHT': 'MF_EIGHT',
    }

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
        module_name=module,
        known_constants=known_constants,
        la_constants_map=la_constants_map,
        la_constants_suffix='_MF',
        module_type_names=frozenset([real_type, complex_type]),
        module_constant_names=frozenset([
            'mf_zero', 'mf_one', 'mf_two', 'mf_half', 'mf_eight',
            'mf_safmin', 'mf_safmax', 'mf_tsml', 'mf_tbig', 'mf_ssml', 'mf_sbig',
            'mf_rtmin', 'mf_rtmax',
            'mf_to_double', 'mf_real',
        ]),
        module_generic_names=frozenset([
            'abs', 'sqrt', 'sin', 'cos', 'tan', 'exp', 'log', 'log10',
            'atan', 'asin', 'acos', 'aint', 'anint', 'sinh', 'cosh', 'tanh',
            'asinh', 'acosh', 'atanh', 'erf', 'erfc', 'erfc_scaled',
            'gamma', 'log_gamma',
            'bessel_j0', 'bessel_j1', 'bessel_y0', 'bessel_y1',
            'bessel_jn', 'bessel_yn',
            'fraction', 'rrspacing', 'spacing', 'epsilon', 'huge', 'tiny',
            'aimag', 'conjg',
            'sign', 'mod', 'atan2', 'dim', 'modulo', 'hypot', 'nearest',
            'min', 'max',
            'maxval', 'minval', 'maxloc', 'minloc', 'sum', 'product',
            'dot_product', 'norm2', 'findloc', 'matmul',
            'digits', 'maxexponent', 'minexponent', 'radix',
            'precision', 'range', 'exponent',
            'dble', 'int', 'nint', 'real', 'cmplx', 'ceiling', 'floor',
            'scale', 'set_exponent',
            'storage_size', 'random_number',
        ]),
        module_operator_generics=(
            'operator(+)', 'operator(-)', 'operator(*)', 'operator(/)',
            'operator(**)',
            'operator(==)', 'operator(/=)',
            'operator(<)', 'operator(>)', 'operator(<=)', 'operator(>=)',
            'assignment(=)',
        ),
        c_real_type='float64x2_t',
        c_complex_type='complex128x2_t',
        c_mpi_real='MPI_FLOAT64X2',
        c_mpi_complex='MPI_COMPLEX128X2',
        c_mpi_sum_real='MPI_DD_SUM',
        c_mpi_sum_complex='MPI_ZZ_SUM',
        c_header='multifloats_bridge.h',
    )
