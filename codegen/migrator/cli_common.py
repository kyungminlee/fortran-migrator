"""Shared CLI helpers used across the migrator subcommands.

Small argument-to-object adapters (target-mode construction, parser
selection) and path/name conventions that several ``cmd_*`` handlers and
``staging`` depend on. Kept in their own module so both ``__main__`` and
``staging`` can import them without a circular dependency. Extracted
verbatim from ``__main__``.
"""
from pathlib import Path

from .target_mode import load_target

# Every non-empty ``la_constants_suffix`` across codegen/targets/*.yaml
# (kind10 → _ey, kind16 → _qx, multifloats → _mw). This is the universe
# of per-target LA_CONSTANTS / LA_XISNAN module suffixes used to route
# the helper pair between OWN/FOREIGN sets during builds and staging.
# The staged top-level CMakeLists.txt (cmake/CMakeLists.txt, LAPACK
# EXCLUDE_REGEX) hardcodes the same universe — keep them in sync.
LA_SUFFIXES = ('_ey', '_qx', '_mw')

# Target used when ``--target`` is absent. Both the argparse default and
# the two ``getattr(args, 'target', None) or ...`` fallbacks (here and in
# ``staging.cmd_stage``, which must name the target before it can build a
# TargetMode) read it from here.
DEFAULT_TARGET = 'kind16'


def dedupe_by_inode(paths) -> list[Path]:
    """Sorted, de-duplicated view of ``paths``, keyed by ``(st_dev, st_ino)``.

    Identity is the inode rather than the name because on a
    case-insensitive filesystem ``*.f`` and ``*.F`` glob the same
    physical file, and a recipe's ``extra_migrate_files`` can name a
    source its ``source_dir`` glob already found. Paths that cannot be
    ``stat``-ed fall back to their string spelling, so a broken symlink
    is kept once rather than dropped.
    """
    seen: dict[tuple, Path] = {}
    for f in paths:
        try:
            st = f.stat()
            key = (st.st_dev, st.st_ino)
        except OSError:
            key = ('missing', str(f))
        seen.setdefault(key, f)
    return sorted(seen.values())


def la_helper_pairs(target_mode) -> tuple[set[str], set[str]]:
    """Return ``(own, foreign)`` LA_CONSTANTS / LA_XISNAN stem sets.

    ``own`` holds the active target's pair-suffixed helper stems;
    ``foreign`` holds every other target's, i.e. the stems the current
    build must not compile from the shared LAPACK SRC directory.
    """
    own_suffix = target_mode.la_constants_suffix.lower()
    own = {f'la_constants{own_suffix}', f'la_xisnan{own_suffix}'}
    foreign = {
        f'la_{base}{sfx}'
        for base in ('constants', 'xisnan')
        for sfx in LA_SUFFIXES
        if sfx != own_suffix
    }
    return own, foreign


def route_sources(files, target_mode, config, independent, *, rel,
                  precision_stems=frozenset(),
                  strip_trailing_underscore: bool = False,
                  copy_files_to: str = 'precision',
                  ) -> tuple[list[str], list[str]]:
    """Split migrated sources into the common and precision archives.

    Returns ``(common, precision)`` in the order ``files`` came in.
    Sources belonging to *another* target's LA_CONSTANTS / LA_XISNAN
    pair are dropped entirely: all three pairs sit in the shared LAPACK
    SRC dir and only the target's own pair belongs in this tree — the
    ``*_mw`` pair does ``use multifloats``, which does not exist in a
    kind10/kind16 build, so compiling it would break the build.

    The routing chain, highest priority first:

    1. ``config.force_common`` — the explicit recipe override that pins
       a stem to the family-independent archive no matter what the
       symbol scanner decided. Rescues files the scanner mis-assigned
       to a non-target family, e.g. the integer BLACS entry points and
       the type-agnostic driver, which are family-independent (one copy
       serves e/y, q/x, m/w) and must not leak into the prefixed
       ``eyblacs`` archive. See ``codegen/recipes/blacs.yaml``.
    2. ``precision_stems`` — stems the caller knows are precision-
       specific by construction (see the parameter, below).
    3. ``config.copy_files``, when ``copy_files_to='precision'``.
    4. ``independent`` — the symbol scanner's verdict.
    5. ``config.copy_files``, when ``copy_files_to='common'``.
    6. everything else → precision.

    Callers differ in four ways, each one a parameter:

    ``rel``
        Formats a :class:`Path` into the manifest spelling. ``cmd_stage``
        writes ``src/<name>``; ``cmd_build`` writes the path relative to
        its output dir. Caller-supplied so neither manifest changes.
    ``precision_stems``
        Upper-cased stems forced to PRECISION ahead of ``independent``.
        ``cmd_stage`` passes the target's *own* LA_CONSTANTS/LA_XISNAN
        pair (single-precision: only this target's E/Y or Q/X block, so
        it belongs in libelapack / libqlapack, never lapack_common)
        together with its precision-rename targets — files the rename
        map gave a precision prefix (``disnan.f`` → ``qisnan.f``) whose
        renamed stem can coincide with a name the classifier marked
        ``independent``, since the split ``la_xisnan_{ey,qx}.F90``
        modules define ``EISNAN``/``QISNAN``/``ELAISNAN``/``QLAISNAN``
        and Q/E/X/Y are not S/D/C/Z. Without the guard the migrated
        plain-F77 isnan externals land in the shared archive and then
        collide across targets in a combined release tree — same
        filename, different content. ``cmd_build`` passes nothing.
    ``strip_trailing_underscore``
        Also test ``force_common`` against the stem with a trailing
        Fortran underscore removed, so a recipe can list the logical
        name (``igamx2d``) exactly as ``skip_files`` does — the C
        migrator strips it, see ``c_migrator.py``, and this keeps the
        two levers' conventions aligned. ``cmd_stage`` wants it (78
        ``*_.c`` stems in BLACS/SRC); ``cmd_build`` must not have it,
        or its ``force_common`` entries would re-route.
    ``copy_files_to``
        Where a surviving ``copy_files`` entry goes. ``cmd_stage``
        sends it to COMMON: the target-specific verbatim copies were
        already claimed by ``precision_stems``, so what is left is
        genuinely shareable (MUMPS common modules, PBLAS Fortran
        helpers) and the symbol scanner may never have visited it — a
        Fortran ``copy_files`` entry in a C recipe never appears in
        ``independent``, so it needs an explicit route. ``cmd_build``
        sends it to PRECISION: that single-target build has no
        per-target special-casing, and some entries USE precision
        modules (LAPACK's ``LA_XISNAN_QX`` USEs ``LA_CONSTANTS_QX``),
        so a common-lib copy would create a forbidden common →
        precision module dependency.
    """
    if copy_files_to not in ('common', 'precision'):
        raise ValueError(
            f"copy_files_to must be 'common' or 'precision', "
            f'got {copy_files_to!r}')

    # Only ``foreign`` is used here. A caller that also wants the
    # target's OWN pair routed to PRECISION passes it back in through
    # ``precision_stems`` — which is what makes the own-pair route
    # opt-in rather than something ``cmd_build`` inherits silently.
    _, foreign = la_helper_pairs(target_mode)
    foreign_upper = {s.upper() for s in foreign}

    common_files: list[str] = []
    precision_files: list[str] = []
    for f in files:
        stem = f.stem.upper()
        if stem in foreign_upper:
            continue
        name = rel(f)
        stem_nodeco = (stem[:-1] if strip_trailing_underscore
                       and stem.endswith('_') else stem)
        if stem in config.force_common or stem_nodeco in config.force_common:
            common_files.append(name)
        elif stem in precision_stems:
            precision_files.append(name)
        elif copy_files_to == 'precision' and stem in config.copy_files:
            precision_files.append(name)
        elif stem in independent:
            common_files.append(name)
        elif copy_files_to == 'common' and stem in config.copy_files:
            common_files.append(name)
        else:
            precision_files.append(name)
    return common_files, precision_files


def recipe_project_root(recipe_path: Path) -> Path:
    """Default project root inferred from a recipe path.

    Recipes live at ``<root>/codegen/recipes/<name>.yaml``, three
    levels below the repo root. The path is used as given — callers
    that want symlink resolution resolve before calling.
    """
    return recipe_path.parent.parent.parent


def migrator_project_root() -> Path:
    """Project root located from this package's position on disk
    (``<root>/codegen/migrator/`` → repo root)."""
    return Path(__file__).resolve().parent.parent.parent


def target_name(args) -> str:
    """Target name from CLI arguments, defaulted.

    Separate from :func:`get_target_mode` because ``cmd_stage`` must
    compare the name against ``BASELINE_TARGETS`` before it knows
    whether a :class:`TargetMode` is wanted at all.
    """
    return getattr(args, 'target', None) or DEFAULT_TARGET


def get_target_mode(args):
    """Construct TargetMode based on CLI arguments."""
    return load_target(target_name(args))


def parser_args(args):
    """Extract parser/parser_cmd from CLI args."""
    parser = getattr(args, 'parser', None)
    parser_cmd = getattr(args, 'parser_cmd', None)
    return parser, parser_cmd
