"""``stage`` subcommand: materialize a self-contained migrated build tree.

Enumerates and migrates each library's sources into a staging directory,
patches the LibSeq MPI shim, classifies symbols, and lays down the CMake
manifest a downstream build consumes. Extracted verbatim from ``__main__``.
"""
import re
import shutil
import sys
from pathlib import Path

from .cli_common import (
    dedupe_by_inode, get_target_mode, la_helper_pairs, migrator_project_root,
    parser_args, target_name,
)
from .libseq_patch import patch_libseq_mpi_f
from .pipeline import classify_recipe_symbols, run_migration
from .prepare import prepare_recipe, run_prepare


# BLACS-style dual-entry-point detector used by ``cmd_stage`` to
# identify C sources that switch their public symbol via the
# ``INTFACE == C_CALL`` / ``CallFromC`` macros. Hoisted to module scope
# so it's compiled once instead of per-library in the cmd_stage loop.
_DUAL_ENTRY_C_RE = re.compile(
    r'#\s*if\s*\(?\s*INTFACE\s*==\s*C_CALL\b'
    r'|#\s*ifdef\s+CallFromC\b'
    r'|#\s*if\s+defined\s*\(\s*CallFromC\s*\)',
)


# --- Shared C struct ABI guard -------------------------------------------
#
# The migrator produces two copies of every shared C header: a verbatim
# type-generic reference under ``_X_src/`` (never migrated) and a typed
# copy under ``X/src/``. Both are linked into the same address space as MKL
# in an MKL-first flatten build, so every struct they define MUST keep a
# byte-identical layout. A single regression already happened: the
# case-insensitive rename conflated the PBTYP_T field ``Cgesd2d`` with the
# BLACS routine ``cgesd2d_``, and the header-duplication pass injected a dead
# ``Wgesd2d`` twin into the typed struct *definition*, shifting the following
# function-pointer slots and NULL-dispatching at runtime (SIGSEGV in
# ``PB_CpaxpbyNN`` at np>=2). Link time cannot catch this — the symbol names
# still match. So we assert it structurally right after staging: every
# aggregate defined in a reference header must survive member-for-member in
# the typed copy. New typedefs the KIND path legitimately adds (QCOMPLEX) are
# ignored — the check only flags reference structs that changed or vanished.

def _decomment_c(text: str) -> str:
    """Blank out C comments and string/char literals so brace-matching and
    struct-body comparison are not fooled by braces or ``struct`` keywords
    inside comments or literals."""
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == '/' and i + 1 < n and text[i + 1] == '*':
            j = text.find('*/', i + 2)
            i = n if j == -1 else j + 2
            out.append(' ')
            continue
        if c == '/' and i + 1 < n and text[i + 1] == '/':
            j = text.find('\n', i)
            i = n if j == -1 else j
            continue
        if c in '"\'':
            q = c
            i += 1
            while i < n:
                if text[i] == '\\':
                    i += 2
                    continue
                if text[i] == q:
                    i += 1
                    break
                i += 1
            out.append(' ')
            continue
        out.append(c)
        i += 1
    return ''.join(out)


def _extract_c_aggregates(text: str) -> dict:
    """Map ``struct|union|enum <name>`` -> whitespace-normalized member body
    for every *named* aggregate definition (one carrying a ``{...}`` body) in
    a C header. Anonymous aggregates are skipped: they cannot be matched
    across the two copies (their only key would be a source offset, which the
    typed copy's inserted bridge-include shifts) and they are not part of any
    cross-TU ABI anyway."""
    code = _decomment_c(text)
    aggregates = {}
    for m in re.finditer(
        r'\b(struct|union|enum)\b[ \t]*([A-Za-z_]\w*)?[ \t\n]*\{', code
    ):
        kind = m.group(1)
        tag = m.group(2) or ''
        brace = m.end() - 1
        depth = 0
        k = brace
        while k < len(code):
            if code[k] == '{':
                depth += 1
            elif code[k] == '}':
                depth -= 1
                if depth == 0:
                    break
            k += 1
        body = code[brace + 1:k]
        semi = code.find(';', k)
        trailer = code[k + 1:semi] if semi != -1 else ''
        names = ([tag] if tag else []) + re.findall(r'[A-Za-z_]\w*', trailer)
        if not names:
            continue  # anonymous aggregate — cannot match reliably, skip
        norm = re.sub(r'\s+', ' ', body).strip()
        for name in names:
            aggregates[f'{kind} {name}'] = norm
    return aggregates


def _assert_shared_struct_abi(staging_dir: Path) -> None:
    """Post-stage guard: every C aggregate defined in a reference ``_X_src``
    header must appear member-for-member identical in the typed ``X/src``
    copy. A divergence means the case-insensitive rename rewrote a shared
    struct layout (the PBTYP_T W-field class of bug) — which links cleanly but
    crashes at runtime. Abort staging loudly instead."""
    violations = []
    for ref_dir in sorted(staging_dir.glob('_*_src')):
        typ_dir = staging_dir / ref_dir.name[1:].removesuffix('_src') / 'src'
        if not typ_dir.is_dir():
            continue
        for ref_h in sorted(ref_dir.glob('*.h')):
            typ_h = typ_dir / ref_h.name
            if not typ_h.is_file():
                continue
            ref_aggs = _extract_c_aggregates(ref_h.read_text(errors='replace'))
            typ_aggs = _extract_c_aggregates(typ_h.read_text(errors='replace'))
            for name, ref_body in ref_aggs.items():
                typ_body = typ_aggs.get(name)
                if typ_body is None:
                    violations.append(
                        f'  {typ_h}: `{name}` present in reference '
                        f'{ref_dir.name}/{ref_h.name} but missing/renamed in '
                        f'typed copy')
                elif typ_body != ref_body:
                    violations.append(
                        f'  {typ_h}: `{name}` layout diverged from '
                        f'{ref_dir.name}/{ref_h.name}\n'
                        f'      reference: {ref_body}\n'
                        f'      typed:     {typ_body}')
    if violations:
        raise SystemExit(
            'staging aborted: shared C struct ABI diverged between the '
            'type-generic reference headers and the migrated typed copies '
            '(PBTYP_T W-field class of bug).\n'
            'A case-insensitive rename likely rewrote a struct member whose '
            'Titlecase spelling collides with a routine name; the typed '
            'archive would link cleanly against MKL/reference plumbing but '
            'NULL-dispatch at runtime. Protect the member in c_migrator.py '
            '(see _PBTYP_C_BINDING_FIELDS).\n\n' + '\n'.join(violations))


# Migrated-source glob patterns, per recipe language. A C recipe's
# staged tree also carries Fortran: ``copy_files`` plants
# precision-independent helpers (PBLAS/SRC/pilaenv.f) that CMake compiles
# into the same archive. A Fortran recipe never emits ``.c``.
_FORTRAN_SOURCE_GLOBS = ('*.f', '*.F', '*.f90', '*.F90')
_SOURCE_GLOBS: dict[str, tuple[str, ...]] = {
    'c': ('*.c', *_FORTRAN_SOURCE_GLOBS),
}


def _collect_source_files(src_dir: Path, language: str) -> list[Path]:
    """Discover migrated source files in ``src_dir`` for the given language.

    Honors all four Fortran extension cases (``.f``/``.F``/``.f90``/``.F90``).
    Dedupe uses ``(st_dev, st_ino)`` so case-insensitive filesystems (where
    ``*.f`` and ``*.F`` glob the same physical file) do not double-stage.
    """
    patterns = _SOURCE_GLOBS.get(language, _FORTRAN_SOURCE_GLOBS)
    return dedupe_by_inode(f for pat in patterns for f in src_dir.glob(pat))


# Topologically sorted library build order for the unified CMake project.
# Each entry is (library_name, recipe_filename).
LIBRARY_ORDER = [
    ('blas',        'blas.yaml'),
    ('xblas',       'xblas.yaml'),
    ('blacs',       'blacs.yaml'),
    ('lapack',      'lapack.yaml'),
    ('ptzblas',     'ptzblas.yaml'),
    # Leaf symbol-provider for NUMROC / ICEIL / ILCM: kept so pbblas.yaml
    # keeps its cycle-free numroc symbol context. Its migrated output is
    # no longer built — the three helpers are folded into scalapack_common
    # (migrated by scalapack.yaml). Still staged (harmless, unused manifest).
    ('scalapack_tools', 'scalapack_tools.yaml'),
    ('pbblas',      'pbblas.yaml'),
    ('pblas',       'pblas.yaml'),
    ('scalapack',   'scalapack.yaml'),
    ('scalapack_c', 'scalapack_c.yaml'),
    ('mumps',       'mumps.yaml'),
]


BASELINE_TARGETS = ('kind4', 'kind8')


# Top-level CMakeLists + its include() modules, FortranCompiler module,
# and the extended-precision detector. Copied verbatim into every staged
# tree (migrated and baseline alike).
_CMAKE_GLUE_FILES = (
    'CMakeLists.txt',
    'EplinalgLibraryHelpers.cmake',
    'EplinalgMumps.cmake',
    'EplinalgInstall.cmake',
    'EplinalgPrivatizeGate.cmake',
    'EplinalgPrivatizeGateCheck.cmake',
    'FortranCompiler.cmake',
    'DetectExtendedPrecision.cmake',
)

# Repo-root files that ride along: ``CMakePresets.json`` so users can
# `cmake --preset=linux-impi` from the staged tree without having to
# re-discover Intel MPI's wrapper paths, and ``VERSION`` so the staged
# CMakeLists.txt's project() call carries the release version.
_ROOT_FILES = ('CMakePresets.json', 'VERSION')


# Standard-precision source directories staged for the std archives
# built alongside each migrated extension. The CMakeLists.txt invokes
# add_standard_fortran_library / add_standard_c_library against these
# directories. Sibling to _refblas_src/_reflapack_src but used for
# production link deps, not just tests.
# Each entry is (dst_name, extern-relative source, recipe name or None).
# The recipe column is consumed only by baseline staging: libraries with
# a recipe stage from build/staged-sources/<lib>/ so the baseline column
# links against the patched archives; libraries without one (TOOLS /
# REDIST / libseq / MUMPS include/) come straight from extern/.
# _pblas_src/ includes PTOOLS/ as a child subdirectory (matching
# the upstream layout) so PTOOLS sources' ``#include "../pblas.h"``
# resolves without an explicit include-path remap. Same shape for
# _scalapack_src/ which contains REDIST/SRC and shares the
# ``../some_header.h`` convention. PBLAS's internal subdirectories
# (PBBLAS, PTZBLAS) are NOT included under _pblas_src — those are
# owned by the separate ptzblas / pbblas std archives.
# _scalapack_src/ is the upstream SRC/ tree; scalapack_c REDIST
# routines live alongside SRC under REDIST/SRC, but we stage them
# together inside _scalapack_src/REDIST/ so REDIST sources'
# ``#include "../redist.h"`` (or the matching SRC headers)
# resolve relative to _scalapack_src/.
_STD_DIRS: list[tuple[str, str, str | None]] = [
    ('_blacs_src',            'scalapack-2.2.3/BLACS/SRC',         'blacs'),
    ('_pblas_src',            'scalapack-2.2.3/PBLAS/SRC',         'pblas'),
    ('_ptzblas_src',          'scalapack-2.2.3/PBLAS/SRC/PTZBLAS', 'ptzblas'),
    ('_pbblas_src',           'scalapack-2.2.3/PBLAS/SRC/PBBLAS',  'pbblas'),
    ('_scalapack_src',        'scalapack-2.2.3/SRC',               'scalapack'),
    ('_scalapack_tools_src',  'scalapack-2.2.3/TOOLS',             None),
    ('_scalapack_redist_src', 'scalapack-2.2.3/REDIST/SRC',        None),
    # MUMPS sequential MPI stub. Lets cmake build a single-process
    # ``libmpiseq`` archive alongside the migrated qxmumps; tests can
    # link it instead of MPI::MPI_Fortran for plain (no mpiexec)
    # executables. Stubs print a "should not be called" error if a
    # collective/comm primitive that requires multi-rank coordination
    # is invoked, so libseq is NPROCS=1-only by construction.
    ('_mpiseq_src',           'MUMPS_5.9.0/libseq',                None),
    # MUMPS upstream src/ + include/. The recipe (which is fortran-
    # only) skips every *MUMPS_C / MUMPS_C_TYPES header and every
    # *.c file, so the migrated qxmumps archive ships without a C
    # interface. tests/mumps's C-bridge build re-uses upstream
    # mumps_c.c (compiled twice with quad-precision type overrides
    # supplied from tests/mumps/c/include/, see B2 in
    # tests/mumps/TODO.md), plus mumps_common.c, mumps_addr.c, and
    # the IO/save/thread/metis/pord/scotch helpers, all of which are
    # type-agnostic and compile verbatim. Staging the whole src/
    # tree (including the .F siblings we don't need here) is
    # cheaper than per-file plumbing and matches the convention.
    ('_mumps_upstream_src',     'MUMPS_5.9.0/src',                 'mumps'),
    ('_mumps_upstream_include', 'MUMPS_5.9.0/include',             None),
    # PORD nested-dissection ordering — ships in-tree with MUMPS and
    # is self-contained standard C (no MPI / external dep). Staging
    # its algorithm sources (PORD/lib) + headers (PORD/include) lets
    # cmake build ``libpord`` and define ``-Dpord`` so ICNTL(7)=4
    # works; without it mumps_pord.c compiles as an inert stub.
    # Precision-agnostic (permutes the integer adjacency graph), so a
    # single build serves every migrated arithmetic.
    ('_mumps_pord_src',         'MUMPS_5.9.0/PORD/lib',            None),
    ('_mumps_pord_include',     'MUMPS_5.9.0/PORD/include',        None),
    # METIS 5.1.0 nested-dissection / k-way ordering — vendored under
    # extern/metis-5.1.0 with every public API symbol privately
    # namespaced (METIS_<X> → METIS_MUMPS_<X>, internal libmetis__ →
    # libmetis_MUMPS_) so this copy can never clash with a system
    # METIS at link time; the MUMPS caller sites were renamed to
    # match. Staging GKlib + libmetis sources and the public header
    # lets cmake build ``libmetis`` and define ``-Dmetis`` so
    # ICNTL(7)=5 works; without it the mumps_metis*.c compile as inert
    # stubs. Integer-graph only, so a single 32-bit-idx build serves
    # every migrated arithmetic.
    ('_mumps_metis_gklib',      'metis-5.1.0/GKlib',               None),
    ('_mumps_metis_lib',        'metis-5.1.0/libmetis',            None),
    ('_mumps_metis_include',    'metis-5.1.0/include',             None),
    # Scotch 7.0.4 sequential ordering (ICNTL(7)=3) — vendored under
    # extern/scotch-7.0.4, built with -DSCOTCH_NAME_SUFFIX=_mumps so
    # every public SCOTCH_* and internal _SCOTCH* symbol carries a
    # _mumps suffix and this copy can never clash with a system Scotch
    # at link time; the MUMPS caller sites resolve to the suffixed
    # names through scotch_rename_mumps.h. The bison/flex parser and
    # scotch.h/scotchf.h are pre-generated and vendored, so the build
    # needs no bison/flex. Staging libscotch + esmumps sources and the
    # generated headers lets cmake build ``scotch``/``esmumps``
    # and define ``-Dscotch`` so ICNTL(7)=3 works; without
    # it the mumps_scotch*.c compile as inert stubs. Integer-graph
    # only, so a single build serves every migrated arithmetic.
    ('_mumps_scotch_libsrc',    'scotch-7.0.4/libscotch',          None),
    ('_mumps_scotch_esmumps',   'scotch-7.0.4/esmumps',            None),
    ('_mumps_scotch_include',   'scotch-7.0.4/include',            None),
]


def _copy_cmake_glue(proj_root: Path, staging_dir: Path,
                     warn: bool = False) -> None:
    """Copy the CMake glue files into the staged tree."""
    cmake_dir = proj_root / 'cmake'
    for cmake_file in _CMAKE_GLUE_FILES:
        src = cmake_dir / cmake_file
        if src.exists():
            shutil.copy2(src, staging_dir / cmake_file)
        elif warn:
            print(f'Warning: {src} not found')


def _copy_root_files(proj_root: Path, staging_dir: Path,
                     warn: bool = False) -> None:
    """Copy the ride-along repo-root files into the staged tree."""
    for root_file in _ROOT_FILES:
        src = proj_root / root_file
        if src.exists():
            shutil.copy2(src, staging_dir / root_file)
        elif warn:
            print(f'Warning: {src} not found')


def _copy_rename_script(proj_root: Path, staging_dir: Path) -> None:
    """Plant the refblas_quad / reflapack_quad symbol-rename helper
    alongside the other build-time scripts so tests/blas/refblas and
    tests/lapack/reflapack can locate it via find_file in the staging
    tree (see those CMakeLists for the search-path list)."""
    scripts_src = proj_root / 'scripts' / 'refquad_rename_archive.sh'
    if scripts_src.exists():
        scripts_dst = staging_dir / 'scripts'
        scripts_dst.mkdir(exist_ok=True)
        shutil.copy2(scripts_src, scripts_dst / scripts_src.name)


def _replace_tree(src: Path, dst: Path, ignore=None) -> bool:
    """Copy the directory ``src`` onto ``dst``, replacing it wholesale.

    Returns False (having done nothing) when ``src`` is not a directory:
    every staging source is optional — a missing upstream tree means the
    build falls back to a system library or simply omits the component.
    """
    if not src.is_dir():
        return False
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=ignore)
    return True


def _copy_tests(proj_root: Path, staging_dir: Path) -> None:
    """Copy tests/ subtree so the unified CMakeLists.txt can pick it up
    via add_subdirectory(tests) when BUILD_TESTING=ON."""
    _replace_tree(proj_root / 'test' / 'integration', staging_dir / 'tests')


def _copy_runtime_mpiseq(proj_root: Path, staging_dir: Path) -> None:
    """Plant runtime/mpiseq (libmpiseq CMakeLists + first-party stub
    sources) into the staged tree, so the parent build's
    add_subdirectory(runtime/mpiseq) resolves. The upstream libseq
    sources it compiles land separately in _mpiseq_src/."""
    _replace_tree(proj_root / 'src' / 'mpiseq', staging_dir / 'runtime' / 'mpiseq')


def _select_libraries(requested) -> list[tuple[str, str]]:
    """The ``(name, recipe_file)`` pairs to stage, in :data:`LIBRARY_ORDER`.

    ``requested`` is the ``--libraries`` filter; empty means all of them.
    An unknown name is a hard error rather than a silent no-op — it is
    almost always a typo, and staging would otherwise look successful
    while producing nothing.
    """
    if not requested:
        return list(LIBRARY_ORDER)
    lib_set = set(requested)
    valid = {n for n, _ in LIBRARY_ORDER}
    unknown = lib_set - valid
    if unknown:
        sys.exit(
            f'error: unknown library name(s) in --libraries: '
            f'{sorted(unknown)}. Valid: {sorted(valid)}'
        )
    return [(n, r) for n, r in LIBRARY_ORDER if n in lib_set]


def _route_manifest_sources(files: list[Path], config, target_mode,
                            classification) -> tuple[list[str], list[str]]:
    """Split staged sources into the common and precision archives.

    Returns ``(common, precision)`` as ``src/<name>`` strings, in the
    order ``files`` came in. Sources belonging to another target's
    LA_CONSTANTS / LA_XISNAN pair are dropped entirely.
    """
    independent = classification.independent

    # Files that the precision-rename map RENAMES are precision-
    # specific by construction — they were given a precision prefix
    # (e.g. ``disnan.f`` → ``qisnan.f``, family ``#RISNAN``, d→q).
    # Their renamed stem can COINCIDE with a name the classifier
    # marked ``independent``: the split ``la_xisnan_ey.F90`` /
    # ``la_xisnan_qx.F90`` modules define procedures ``EISNAN`` /
    # ``QISNAN`` / ``ELAISNAN`` / ``QLAISNAN``, which land in
    # ``independent`` because Q/E/X/Y are not S/D/C/Z.
    # Without this guard the migrated plain-F77 isnan externals
    # (``qisnan_``/``eisnan_`` — whose bodies differ per target) get
    # mis-filed into the shared ``lapack_common`` archive and then
    # collide across targets in a combined release tree (same
    # filename, different content). A rename target is never a
    # legitimate common member, so route it to PRECISION regardless
    # of the polluted independent set.
    rename_targets = {
        t.upper()
        for t in classification.build_rename_map(target_mode).values()
    }

    # Per-target LA_CONSTANTS / LA_XISNAN precision helpers. Each
    # extended target owns exactly one prefix-pair module —
    #   kind10 → *_ey   kind16 → *_qx   multifloats → *_mw
    # (suffix taken from the target's ``la_constants_suffix``). All
    # three pairs sit in the shared LAPACK SRC dir; only the
    # target's own pair belongs in this staging tree — the other
    # targets' pairs are excluded, so e.g. the *_ey/*_qx SRC files
    # stay out of a multifloats build.
    la_own, la_foreign = la_helper_pairs(target_mode)

    common_files, precision_files = [], []
    for f in files:
        if f.stem in la_foreign:
            continue
        rel = f'src/{f.name}'
        stem = f.stem.upper()
        # ``force_common`` is the explicit recipe override that pins a
        # stem to the family-independent ``_common`` archive no matter
        # what the symbol scanner decided. It takes priority over every
        # other route (including the PRECISION defaults below) so a
        # recipe author can rescue files the scanner mis-assigned to a
        # non-target family — e.g. the integer BLACS entry points and
        # the type-agnostic driver, which are family-independent (one
        # copy serves e/y, q/x, m/w) and must not leak into the
        # prefixed ``eyblacs`` archive. See ``codegen/recipes/blacs.yaml``.
        # Match with OR without the trailing Fortran-underscore so a
        # recipe can list the logical name (``igamx2d``) exactly as
        # ``skip_files`` does (the C migrator strips it — see
        # c_migrator.py), keeping the two levers' conventions aligned.
        stem_nodeco = stem[:-1] if stem.endswith('_') else stem
        if stem in config.force_common or stem_nodeco in config.force_common:
            common_files.append(rel)
        # The target's own LA_CONSTANTS/LA_XISNAN pair is single-
        # precision (only this target's E/Y or Q/X block), so it is
        # precision-specific: route it to PRECISION so it lands in
        # the prefixed archive (libelapack / libqlapack), never the
        # shared lapack_common. A precision-rename target (e.g. the
        # migrated qisnan.f, whose renamed stem may collide with an
        # ``independent`` module-procedure name) is likewise
        # precision-specific by construction — see ``rename_targets``.
        elif f.stem in la_own or stem in rename_targets:
            precision_files.append(rel)
        # Any ``copy_files`` entry that survives to this branch is
        # precision-independent (staged verbatim, no prefix rename,
        # shareable across targets): the target-specific verbatim
        # copies — LAPACK's LA_CONSTANTS/LA_XISNAN pairs — were
        # already routed by the ``la_own``/``la_foreign`` checks
        # above. The symbol scanner may never have visited these —
        # e.g. a Fortran ``copy_files`` entry in a C recipe — so
        # they won't appear in ``independent``. Treat them as
        # common explicitly. (``cmd_build`` in __main__.py lacks
        # the LA_* special-casing and therefore routes copy_files
        # to PRECISION instead — see the comment there.)
        elif stem in independent or stem in config.copy_files:
            common_files.append(rel)
        else:
            precision_files.append(rel)
    return common_files, precision_files


def _detect_dual_entry_sources(files: list[Path], config) -> list[str]:
    """C sources that gate their entry-point signature on the
    ``INTFACE == C_CALL`` macro (upstream BLACS pattern).

    Each such file exposes a Fortran-callable symbol (e.g.
    ``blacs_gridinfo_``) in the default build and a C-callable symbol
    (``Cblacs_gridinfo``) when compiled with ``-DCallFromC``. The CMake
    helper compiles these sources twice so the final static library
    ships both entry points. Detection is a cheap regex scan of the
    staged source; any library that does not use the pattern emits an
    empty list.
    """
    if config.language != 'c':
        return []
    dual_files = []
    for f in files:
        if f.suffix.lower() != '.c':
            continue
        try:
            text = f.read_text(errors='replace')
        except OSError:
            continue
        if _DUAL_ENTRY_C_RE.search(text):
            dual_files.append(f'src/{f.name}')
    return dual_files


def _write_manifest(lib_dir: Path, lib_name: str, language: str,
                    common_files: list[str], precision_files: list[str],
                    dual_files: list[str]) -> None:
    """Write a library's ``manifest.cmake`` — the source lists the staged
    tree's ``add_migrated_*`` helpers read."""
    common_list = '\n    '.join(common_files) if common_files else ''
    precision_list = '\n    '.join(precision_files) if precision_files else ''
    dual_list = '\n    '.join(dual_files) if dual_files else ''
    manifest = f"""\
set({lib_name}_COMMON_SOURCES
    {common_list}
)

set({lib_name}_PRECISION_SOURCES
    {precision_list}
)

set({lib_name}_DUAL_INTERFACE_SOURCES
    {dual_list}
)

set({lib_name}_LANGUAGE {language})
"""
    (lib_dir / 'manifest.cmake').write_text(manifest)


def _stage_one_library(lib_name: str, recipe_path: Path, lib_dir: Path,
                       target_mode, proj_root: Path,
                       parser: str | None, parser_cmd: str | None) -> bool:
    """Migrate one library into ``lib_dir`` and write its manifest.

    Returns whether the library went through the ep_ symbol-privatization
    pass, which the caller needs for ``PRIVATIZED_LIBRARIES``.
    """
    src_dir = lib_dir / 'src'
    src_dir.mkdir(parents=True, exist_ok=True)

    print(f'\n{"=" * 60}')
    print(f'  Migrating: {lib_name}')
    print(f'{"=" * 60}')

    run_migration(
        recipe_path=recipe_path,
        output_dir=src_dir,
        target_mode=target_mode,
        dry_run=False,
        project_root=proj_root,
        parser=parser,
        parser_cmd=parser_cmd,
    )

    config = prepare_recipe(recipe_path, proj_root)
    classification = classify_recipe_symbols(config)

    # Pick up precision-independent Fortran helpers staged via
    # ``copy_files`` (e.g. PBLAS/SRC/pilaenv.f) when the recipe is C —
    # CMake's ``add_library(… STATIC …)`` handles mixed C + Fortran
    # sources because both languages are enabled at the top
    # project(). Copy-files are precision-independent by contract,
    # so they land in COMMON_SOURCES.
    files = _collect_source_files(src_dir, config.language)

    common_files, precision_files = _route_manifest_sources(
        files, config, target_mode, classification)
    dual_files = _detect_dual_entry_sources(files, config)

    _write_manifest(lib_dir, lib_name, config.language,
                    common_files, precision_files, dual_files)
    print(f'  Manifest: {len(common_files)} common, '
          f'{len(precision_files)} precision files')
    return bool(config.privatize_symbols)


def _prior_library_list(staging_dir: Path, var_name: str) -> list[str]:
    """Entries of ``var_name`` in an existing ``target_config.cmake``
    whose library directory is still present in the staging tree."""
    prior: list[str] = []
    prior_config = staging_dir / 'target_config.cmake'
    if prior_config.exists():
        m = re.search(
            r'^\s*set\s*\(\s*' + var_name + r'\s+([^)]*)\)',
            prior_config.read_text(),
            re.MULTILINE,
        )
        if m:
            for tok in m.group(1).replace(';', ' ').split():
                tok = tok.strip().strip('"')
                if tok and (staging_dir / tok).is_dir():
                    prior.append(tok)
    return prior


def _merge_staged_lists(staging_dir: Path, staged: list[str],
                        privatized: list[str]) -> tuple[str, str]:
    """The ``STAGED_LIBRARIES`` / ``PRIVATIZED_LIBRARIES`` values to write.

    Re-staging with a subset of --libraries must not shrink the unified
    build's library list, so both lists union this run's result with
    whatever a prior ``target_config.cmake`` recorded and is still on
    disk. A prior *privatized* entry survives only if its library was NOT
    restaged this run — a restage recomputes the flag from the recipe.
    """
    def merged(prior: list[str], fresh: list[str]) -> str:
        out: list[str] = []
        for n in prior + fresh:
            if n not in out:
                out.append(n)
        return ';'.join(out)

    staged_list = merged(_prior_library_list(staging_dir, 'STAGED_LIBRARIES'),
                         staged)
    privatized_list = merged(
        [n for n in _prior_library_list(staging_dir, 'PRIVATIZED_LIBRARIES')
         if n not in staged],
        privatized)
    return staged_list, privatized_list


def _stage_reference_sources(proj_root: Path, staging_dir: Path) -> None:
    """Stage the vendored Netlib BLAS and LAPACK sources the differential
    precision tests build their references from.

    ``_refblas_src`` feeds refblas_quad (compiled with gfortran's
    -freal-8-real-16 to promote KIND=8 entities to KIND=16 in-place);
    tests fall back to system -lblas if it is absent. ``_reflapack_src``
    does the same for tests/lapack/reflapack/, giving a KIND=16 reference
    to compare the migrated qxlapack/eylapack/mwlapack against. LAPACK
    additionally needs {s,d}lamch.f / {s,d}roundup_lwork.f, which its SRC
    routines call but which live in INSTALL/ — copy them in alongside so
    a single glob compiles the full reference. Both the single- and
    double-precision variants are needed: the genuine standard ``lapack``
    archive carries every arithmetic (the genuine single-precision
    solvers smumps/cmumps call slamch_ / sroundup_lwork_, the double ones
    dlamch_ / droundup_lwork_). These verbatim upstream files are
    unpromoted here — only the migrated LAPACK archive elides
    roundup_lwork via _strip_roundup_lwork; the standard one keeps it.
    """
    netlib = proj_root / 'extern' / 'lapack-3.12.1'
    _replace_tree(netlib / 'BLAS' / 'SRC', staging_dir / '_refblas_src')

    reflapack_dst = staging_dir / '_reflapack_src'
    if _replace_tree(netlib / 'SRC', reflapack_dst):
        install_src = netlib / 'INSTALL'
        for fname in ('slamch.f', 'dlamch.f',
                      'sroundup_lwork.f', 'droundup_lwork.f'):
            src = install_src / fname
            if src.is_file():
                shutil.copy2(src, reflapack_dst / fname)


def cmd_stage(args):
    """Migrate all libraries into a structured staging directory.

    Produces a self-contained directory that can be built with:
        cmake -S <staging> -B <staging>/build && cmake --build <staging>/build -j

    For ``--target kind4`` / ``--target kind8`` the migration step is
    skipped: those are un-migrated single/double precision baselines used
    for diff'ing against the quad reference, so all that's needed is a
    staging tree with upstream sources, tests/, and a target_config.cmake
    that points the test framework at the standard archive (LIB_PREFIX="").
    """
    target_str = target_name(args)
    if target_str in BASELINE_TARGETS:
        return _stage_baseline(args, target_str)

    staging_dir = args.staging_dir.resolve()
    target_mode = get_target_mode(args)
    parser, parser_cmd = parser_args(args)
    proj_root = args.project_root or migrator_project_root()
    recipes_dir = proj_root / 'codegen' / 'recipes'

    libraries = _select_libraries(args.libraries)

    # The real / complex family letters (e/y, q/x, m/w). Stamped into
    # target_config.cmake as LIB_PREFIX / LIB_PREFIX_COMPLEX (symbol
    # prefixes) and, concatenated, as LIB_PAIR_PREFIX — the target /
    # archive / package name prefix (eyblas, qxlapack, mwmumps, ...).
    pmap = target_mode.prefix_map
    lib_prefix = pmap['R'].lower()
    lib_prefix_complex = pmap['C'].lower()

    staged = []
    privatized = []
    for lib_name, recipe_file in libraries:
        recipe_path = recipes_dir / recipe_file
        if not recipe_path.exists():
            print(f'Warning: recipe {recipe_path} not found, skipping {lib_name}')
            continue

        privatized_lib = _stage_one_library(
            lib_name, recipe_path, staging_dir / lib_name, target_mode,
            proj_root, parser, parser_cmd)
        staged.append(lib_name)
        if privatized_lib:
            privatized.append(lib_name)

    # Copy MF helper modules into staging so it's self-contained
    needs_mf = target_mode.module_name is not None
    # kind16 (REAL(KIND=16) / __float128) keeps the *standard* MPI datatypes
    # (MPI_REAL16 / MPI_COMPLEX32) but routes every reduction through custom
    # combine ops, because Intel MPI's builtin MPI_SUM/MAX/MIN have no
    # 16-byte-real kernel. Those ops live in a small plain-C helper
    # (quad_mpi.c) plus a Fortran shim (quad_mpi_f.f90) -- the quad analogue
    # of the multifloats bridge, but far simpler (no FetchContent, no C++,
    # no custom datatypes). A KIND-based target with a custom MPI module and
    # no multifloats module_name is exactly kind16.
    needs_quad_mpi = (
        target_mode.c_mpi_module is not None and not needs_mf
    )

    staged_list, privatized_list = _merge_staged_lists(
        staging_dir, staged, privatized)

    if needs_mf:
        # The first-party multifloats-mpi bridge library. It is a
        # standalone library under runtime/ with its own CMakeLists.txt
        # (add_subdirectory'd from cmake/CMakeLists.txt); staging just
        # plants the whole directory so the relative
        # add_subdirectory(runtime/multifloats-mpi) resolves in the staged
        # tree. (The MPICH_SKIP_MPICXX / OMPI_SKIP_MPICXX guard is baked
        # into multifloats_bridge.h itself.)
        _replace_tree(proj_root / 'src' / 'multifloats-mpi',
                      staging_dir / 'runtime' / 'multifloats-mpi')

    if needs_quad_mpi:
        # The first-party quad-mpi bridge library (custom MPI reduce ops
        # on the standard MPI_REAL16 / MPI_COMPLEX32). Like
        # multifloats-mpi it is a standalone library under runtime/ with
        # its own CMakeLists.txt, add_subdirectory'd from the parent build.
        _replace_tree(proj_root / 'src' / 'quad-mpi',
                      staging_dir / 'runtime' / 'quad-mpi')

    # Copy the libmpiseq assembly (CMakeLists + first-party stub sources)
    # wholesale, for every target: the sequential MPI stub library under
    # runtime/mpiseq/ builds against the upstream libseq sources staged
    # into _mpiseq_src/ (add_subdirectory'd from the parent build, gated
    # there on _mpiseq_src presence + Intel MPI).
    _copy_runtime_mpiseq(proj_root, staging_dir)

    _write_target_config(
        staging_dir, target_str=target_mode.name,
        lib_prefix=lib_prefix, lib_prefix_complex=lib_prefix_complex,
        needs_multifloats=needs_mf, needs_quad_mpi=needs_quad_mpi,
        staged_list=staged_list, privatized_list=privatized_list,
        annotate=True)

    # Copy CMake files to staging directory.
    _copy_cmake_glue(proj_root, staging_dir, warn=True)

    # Repo-root ride-along files (presets + VERSION, see _ROOT_FILES).
    _copy_root_files(proj_root, staging_dir, warn=True)

    # The privatize gate audits each PRIVATIZED_LIBRARIES archive against
    # the same checked-in manifest the migrator's rename pass consumed;
    # both copies come from this one staging run, so pass and gate can
    # never disagree within a staged tree.
    if privatized_list:
        manifest_src = proj_root / 'codegen' / 'recipes' / 'privatize_ep_symbols.txt'
        if manifest_src.exists():
            shutil.copy2(manifest_src,
                         staging_dir / 'privatize_ep_symbols.txt')
        else:
            print(f'Warning: {manifest_src} not found — the privatize '
                  'gate will fail at configure time')

    _copy_rename_script(proj_root, staging_dir)

    _copy_tests(proj_root, staging_dir)

    _stage_reference_sources(proj_root, staging_dir)

    # Stage standard-precision source directories for the std archives
    # built alongside each migrated extension (see the _STD_DIRS table
    # for the directory list and its layout rationale). Migrated targets
    # always copy from extern/; the recipe column only matters for
    # baseline staging.
    for dst_name, rel_src, _recipe in _STD_DIRS:
        _replace_tree(proj_root / 'extern' / rel_src, staging_dir / dst_name)

    # libseq's mpi.f bundles BLACS/ScaLAPACK forwarders (which collide
    # with the real migrated archives) and only knows the standard MPI
    # datatypes in MUMPS_COPY (no MPI_REAL16 / MPI_COMPLEX32 needed by
    # qxmumps reductions). Patch the staged copy to fix both. Upstream
    # extern/ stays read-only.
    mpiseq_dst = staging_dir / '_mpiseq_src' / 'mpi.f'
    if mpiseq_dst.is_file():
        patch_libseq_mpi_f(mpiseq_dst)

    # Fail loudly if a case-insensitive rename corrupted any shared C struct
    # layout (raises SystemExit on violation; see the guard's docstring).
    _assert_shared_struct_abi(staging_dir)

    print(f'\n{"=" * 60}')
    print(f'  Staging complete: {len(staged)} libraries')
    print(f'{"=" * 60}')
    print(f'  Target:  {target_mode.name} (prefix: {lib_prefix})')
    print(f'  Output:  {staging_dir}')
    print('\nTo build:')
    print(f'  cmake -S {staging_dir} -B {staging_dir}/build -DCMAKE_BUILD_TYPE=Release')
    print(f'  cmake --build {staging_dir}/build -j')


def _cmake_bool(value) -> str:
    """CMake spelling of a Python truth value."""
    return 'TRUE' if value else 'FALSE'


def _write_target_config(staging_dir: Path, *, target_str: str,
                         subtitle: str | None = None,
                         lib_prefix: str = '', lib_prefix_complex: str = '',
                         needs_multifloats: bool = False,
                         needs_quad_mpi: bool | None = None,
                         staged_list: str = '', privatized_list: str = '',
                         annotate: bool = False) -> None:
    """Write the staged tree's ``target_config.cmake``.

    One writer for both stagings. The defaults are the baseline shape —
    empty prefixes, no multifloats, no staged or privatized archives —
    so ``_stage_baseline`` only has to name the target.

    ``needs_quad_mpi=None`` omits the ``NEEDS_QUAD_MPI`` line entirely:
    the baseline never links quad MPI and its config has never carried
    the variable. ``annotate`` adds the prose comment above
    ``PRIVATIZED_LIBRARIES``, which is likewise only meaningful once
    there are migrated archives to privatize. Both flags exist to keep
    the two files byte-identical to what the separate writers produced.
    """
    lines = [f'# Generated by: python -m migrator stage --target {target_str}']
    if subtitle:
        lines.append(subtitle)
    lines += [
        f'set(TARGET_NAME "{target_str}")',
        f'set(LIB_PREFIX "{lib_prefix}")',
        f'set(LIB_PREFIX_COMPLEX "{lib_prefix_complex}")',
        f'set(LIB_PAIR_PREFIX "{lib_prefix}{lib_prefix_complex}")',
        f'set(NEEDS_MULTIFLOATS {_cmake_bool(needs_multifloats)})',
    ]
    if needs_quad_mpi is not None:
        lines.append(f'set(NEEDS_QUAD_MPI {_cmake_bool(needs_quad_mpi)})')
    # C sources are compiled as C++ exactly when the target needs the
    # multifloats module — the C++ type is what the C code binds to.
    lines += [
        f'set(C_AS_CXX {_cmake_bool(needs_multifloats)})',
        f'set(STAGED_LIBRARIES {staged_list})',
    ]
    if annotate:
        lines += [
            '# Libraries whose migrated sources went through the ep_ symbol-',
            '# privatization pass (task 44): their shared ``_common`` archives'
            ' carry',
            '# the hermetic ep_ engine and install under ep-prefixed filenames'
            ' /',
            '# packages (see EplinalgInstall.cmake).',
        ]
    lines.append(f'set(PRIVATIZED_LIBRARIES {privatized_list})')
    (staging_dir / 'target_config.cmake').write_text(
        '\n'.join(lines) + '\n')


def _stage_baseline(args, target_str: str):
    """Stage an unmigrated baseline tree for kind4 / kind8.

    No per-library migration is run — kind4 / kind8 are the upstream
    S/D/C/Z entry points themselves, served by the ``add_standard_*``
    archives the unified CMake build always assembles. We just need the
    upstream source trees, the tests/ subtree, the cmake/ glue, and a
    target_config.cmake that signals "no migrated archives" via an
    empty LIB_PREFIX so the test framework links the std archive
    directly.
    """
    staging_dir = args.staging_dir.resolve()
    proj_root = args.project_root or migrator_project_root()
    staging_dir.mkdir(parents=True, exist_ok=True)

    print(f'\n{"=" * 60}')
    print(f'  Baseline staging: {target_str} (no migration)')
    print(f'{"=" * 60}')

    # target_config.cmake: empty LIB_PREFIX, no multifloats, no migrated
    # libs in STAGED_LIBRARIES. The parent CMakeLists.txt's
    # add_migrated_* helpers are no-ops without per-lib manifest.cmake
    # files (which we don't write here), so the build resolves to just
    # the standard archives.
    _write_target_config(
        staging_dir, target_str=target_str,
        subtitle=(f'# Baseline (un-migrated) target — see '
                  f'codegen/targets/{target_str}.yaml.'))

    # CMake glue (top-level CMakeLists + its include() modules,
    # FortranCompiler module, presets). Same files cmd_stage copies;
    # baseline reuses them.
    _copy_cmake_glue(proj_root, staging_dir)

    # Repo-root presets + VERSION, same as cmd_stage.
    _copy_root_files(proj_root, staging_dir)

    # libmpiseq assembly (CMakeLists + first-party stub sources), same
    # subtree cmd_stage copies.
    _copy_runtime_mpiseq(proj_root, staging_dir)

    _copy_rename_script(proj_root, staging_dir)

    # tests/ subtree.
    _copy_tests(proj_root, staging_dir)

    # Upstream BLAS / LAPACK / ScaLAPACK / PBLAS / BLACS / MUMPS sources.
    # For libraries with a recipe we stage from build/staged-sources/<lib>/
    # so the baseline column links against the patched archives (closing
    # the gap the kind4/kind8 column previously had — patched in migrated
    # archive but broken in baseline). Libraries without a recipe
    # (TOOLS / REDIST / libseq / MUMPS include/) come straight from
    # extern/.
    recipes_dir = proj_root / 'codegen' / 'recipes'

    def _staged_or_external(rel_src: str, recipe_name: str | None) -> Path:
        if recipe_name:
            recipe_path = recipes_dir / f'{recipe_name}.yaml'
            if recipe_path.exists():
                return run_prepare(recipe_path, project_root=proj_root)
        return proj_root / 'extern' / rel_src

    def _stage_dst(dst_name: str, src: Path) -> None:
        _replace_tree(src, staging_dir / dst_name,
                      ignore=shutil.ignore_patterns('.prepared.stamp'))

    _stage_dst('_refblas_src',
               _staged_or_external('lapack-3.12.1/BLAS/SRC', 'blas'))

    reflapack_src = _staged_or_external('lapack-3.12.1/SRC', 'lapack')
    _stage_dst('_reflapack_src', reflapack_src)
    if reflapack_src.is_dir():
        # Pull dlamch / slamch / droundup_lwork / sroundup_lwork from
        # extern/INSTALL (these aren't in any recipe's source_dir, so
        # they don't get staged; baseline needs them alongside SRC).
        install_src = proj_root / 'extern' / 'lapack-3.12.1' / 'INSTALL'
        reflapack_dst = staging_dir / '_reflapack_src'
        for fname in ('dlamch.f', 'droundup_lwork.f', 'slamch.f',
                      'sroundup_lwork.f'):
            src = install_src / fname
            if src.is_file():
                shutil.copy2(src, reflapack_dst / fname)

    for dst_name, rel_src, recipe_name in _STD_DIRS:
        _stage_dst(dst_name, _staged_or_external(rel_src, recipe_name))

    print(f'  Output:  {staging_dir}')
    print('\nTo build:')
    print(f'  cmake -S {staging_dir} -B {staging_dir}/build -DCMAKE_BUILD_TYPE=Release')
    print(f'  cmake --build {staging_dir}/build -j')
    print(f'  ctest --test-dir {staging_dir}/build')
