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
    # Additional C source directories to *migrate* (not just scan) in
    # the same generic-rename-map pass as ``source_dir``. Used by PBLAS
    # to pull in the PTOOLS/ helper sources alongside the SRC/ entry
    # points so the cloned ddgemm.c entry points have a real
    # PB_Cddtypeset implementation to call. Files are flat-copied into
    # ``output_dir`` (no subdirectory mirroring) so include resolution
    # stays simple.
    extra_c_dirs: list[Path] = field(default_factory=list)
    # Additional Fortran source directories whose files are migrated in
    # the same pass as ``source_dir``. Used by MUMPS to pull in the
    # per-arithmetic header files under ``external/MUMPS_5.8.2/include/``
    # (``dmumps_struc.h`` etc.), which are Fortran content despite the
    # ``.h`` extension and must be migrated so the ``INCLUDE`` statements
    # in ``dmumps_struc_def.F`` resolve against the renamed target file.
    # Files are flat-copied into ``output_dir`` (no subdir mirroring).
    extra_fortran_dirs: list[Path] = field(default_factory=list)
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
    # Target-gated verbatim file overrides. Structure:
    #
    #     overrides:
    #       <target_name>:
    #         src_dir: <path relative to recipe file>
    #         files:
    #           - <filename>
    #           - ...
    #
    # For the active target, each listed file is copied verbatim from
    # ``<recipe_dir>/<src_dir>/<filename>`` to ``<output_dir>/<filename>``
    # after the main C migration has produced clones and header patches,
    # so the override wins. Used for hand-written replacement kernels
    # that cannot be produced by regex substitution (e.g. BI_*vv* for
    # multifloats double-double arithmetic).
    overrides: dict = field(default_factory=dict)
    # Directory containing the recipe file, used to resolve paths in
    # ``overrides`` and similar recipe-relative references.
    recipe_dir: Path | None = None
    # Per-line "keep-kind" manifest: for each source filename (basename,
    # not stem — e.g. ``dana_aux.F``), the set of 1-based line numbers
    # whose ``DOUBLE PRECISION`` declarations must NOT be promoted.
    # Used by MUMPS, where ``DOUBLE PRECISION`` overloads "working
    # precision" and "arithmetic-agnostic DP" (timing, flop counters,
    # MPI_WTIME buffers, Scotch ABI). Generated by
    # ``scripts/mumps_sweep_keep_kind.py``. See ``recipes/README.md``.
    keep_kind_lines: dict[str, frozenset[int]] = field(default_factory=dict)


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

    # Load the keep-kind manifest if specified. The manifest is a plain
    # text file with one ``<path>:<lineno>`` entry per line; the path is
    # ignored (we key by basename) and the lineno is 1-based.
    keep_kind_lines: dict[str, set[int]] = {}
    manifest_rel = data.get('keep_kind_manifest')
    if manifest_rel:
        manifest_path = recipe_path.parent / manifest_rel
        for entry in manifest_path.read_text().splitlines():
            entry = entry.strip()
            if not entry or entry.startswith('#'):
                continue
            path_str, _, lineno_str = entry.rpartition(':')
            if not path_str or not lineno_str:
                continue
            basename = Path(path_str).name
            keep_kind_lines.setdefault(basename, set()).add(int(lineno_str))
    keep_kind_frozen = {k: frozenset(v) for k, v in keep_kind_lines.items()}

    return RecipeConfig(
        library=data['library'],
        language=data['language'],
        source_dir=source_dir,
        extensions=[e.lower() for e in data.get('extensions', ['.f', '.f90'])],
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
        overrides=dict(data.get('overrides') or {}),
        recipe_dir=recipe_path.parent,
        extra_c_dirs=[project_root / d
                      for d in (data.get('extra_c_dirs') or [])],
        extra_fortran_dirs=[project_root / d
                            for d in (data.get('extra_fortran_dirs') or [])],
        keep_kind_lines=keep_kind_frozen,
    )
