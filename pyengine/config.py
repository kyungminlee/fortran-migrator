"""Recipe configuration loader.

Recipes are YAML files that describe a library's source layout and
migration parameters. The engine reads a recipe and performs the
appropriate migration.
"""

from dataclasses import dataclass, field
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore[assignment]


@dataclass
class RecipeConfig:
    """Configuration for a single library migration."""
    library: str
    language: str                     # "fortran" or "c"
    source_dir: Path
    extensions: list[str]
    symbols_method: str = 'scan_source'  # "scan_source" or "nm_library"
    library_path: Path | None = None
    prefix_style: str = 'direct'      # "direct" or "scalapack"
    skip_files: set[str] = field(default_factory=set)
    copy_files: set[str] = field(default_factory=set)  # Copy unchanged (multi-precision utilities)
    # Source stems (uppercase, no extension) whose migrated output
    # should win as canonical, overriding the default D/Z-first
    # preference. Used to route around upstream bugs that live only
    # in the D or Z half of a precision pair (e.g. ScaLAPACK's
    # PZUNGQL / PZUNML2 call PB_TOPGET where they should call
    # PB_TOPSET; PCUNGQL / PCUNML2 have the correct restore).
    prefer_source: set[str] = field(default_factory=set)
    # Local-variable renames applied to the S/C half of each pair
    # *before* light-normalized convergence comparison. Each entry
    # maps an S/C-half identifier to its D/Z-half counterpart (e.g.
    # ``CR: ZR``, ``SX: DX``, ``STEMP: DTEMP``). Local variables are
    # not in the symbol table, so the migrator can't rename them —
    # but they are upstream S/D/C/Z convention drift, not semantic
    # differences. Applied with word boundaries to the in-memory
    # migrated S/C text; the on-disk canonical is never modified.
    local_renames: dict[str, str] = field(default_factory=dict)
    copy_all_originals: bool = False  # For C: copy all files, then add clones
    patches: list[str] = field(default_factory=list)
    depends: list[Path] = field(default_factory=list)  # Dependency recipe paths
    extra_symbol_dirs: list[Path] = field(default_factory=list)  # Extra dirs to scan for symbols
    # Additional C return types to recognize when scanning for function
    # definitions, as regex fragments (e.g. ``r'PBTYP_T\s*\*'``). Used
    # only when ``language == 'c'``; the default set in
    # ``symbol_scanner._C_DEFAULT_RETURN_TYPES`` is always included.
    c_return_types: list[str] = field(default_factory=list)
    # Extra library-specific C typedef renames applied after the
    # built-in double/float/SCOMPLEX/DCOMPLEX substitutions. Each entry
    # has ``from`` (list of source identifiers) and ``to`` (target
    # identifier). Both the ``to`` field and the inserted text in
    # ``header_patches`` support template substitution with the C
    # migrator's template_vars (``{REAL_TYPE}``, ``{COMPLEX_TYPE}``,
    # ``{C_REAL_TYPE}``, ``{RP}``, ``{CP}``, ``{RPU}``, ``{CPU}``).
    c_type_aliases: list[dict] = field(default_factory=list)
    # Insert new content into migrated headers after a literal anchor
    # line. Each entry: ``{'file': <relative path under source_dir>,
    # 'after': <anchor line>, 'insert': <text>}``. Used to define
    # library-specific extended-precision typedefs referenced by
    # c_type_aliases targets.
    header_patches: list[dict] = field(default_factory=list)


def load_recipe(recipe_path: Path,
                project_root: Path | None = None) -> RecipeConfig:
    """Load a recipe YAML file.

    Args:
        recipe_path: Path to the .yaml recipe file.
        project_root: Project root for resolving relative paths.
            Defaults to recipe_path's grandparent.
    """
    if yaml is None:
        raise ImportError(
            'PyYAML is required to load recipe configs. '
            'Install it with: pip install pyyaml'
        )

    if project_root is None:
        # recipes/*.yaml → project root is two levels up
        project_root = recipe_path.parent.parent

    with open(recipe_path) as f:
        data = yaml.safe_load(f)

    source_dir = project_root / data['source_dir']

    library_path = None
    symbols_cfg = data.get('symbols', {})
    if symbols_cfg.get('library_path'):
        library_path = project_root / symbols_cfg['library_path']

    skip = set(s.upper() for s in data.get('skip_files', []))
    copy = set(s.upper() for s in data.get('copy_files', []))
    prefer = set(s.upper() for s in data.get('prefer_source', []))
    local_renames = {
        str(k).upper(): str(v).upper()
        for k, v in (data.get('local_renames') or {}).items()
    }

    # Resolve dependency recipe paths relative to the recipe directory
    depends_raw = data.get('depends', [])
    depends = [recipe_path.parent / d for d in depends_raw]

    # Resolve extra symbol directories relative to the project root
    extra_dirs_raw = data.get('extra_symbol_dirs', [])
    extra_symbol_dirs = [project_root / d for d in extra_dirs_raw]

    return RecipeConfig(
        library=data['library'],
        language=data['language'],
        source_dir=source_dir,
        extensions=data.get('extensions', ['.f', '.f90']),
        symbols_method=symbols_cfg.get('method', 'scan_source'),
        library_path=library_path,
        prefix_style=data.get('prefix', {}).get('style', 'direct'),
        skip_files=skip,
        copy_files=copy,
        prefer_source=prefer,
        local_renames=local_renames,
        copy_all_originals=data.get('copy_all_originals', False),
        patches=data.get('patches', []),
        depends=depends,
        extra_symbol_dirs=extra_symbol_dirs,
        c_return_types=list(data.get('c_return_types', [])),
        c_type_aliases=list(data.get('c_type_aliases', [])),
        header_patches=list(data.get('header_patches', [])),
    )
