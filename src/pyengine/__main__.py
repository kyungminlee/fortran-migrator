"""CLI entry point for the general-purpose migration pipeline.

Usage (from the src/ directory):
    uv run python -m pyengine migrate  ../recipes/blas.yaml output/ --target kind16
    uv run python -m pyengine verify   output/
    uv run python -m pyengine build    ../recipes/blas.yaml output/ --target kind16
    uv run python -m pyengine run      ../recipes/blas.yaml work/ --target kind16
    uv run python -m pyengine stage    /tmp/staging --target multifloats
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from .config import load_recipe
from .pipeline import (
    run_convergence_report, run_divergence_report, run_migration,
)
from .prefix_classifier import classify_symbols
from .symbol_scanner import scan_symbols
from .target_mode import load_target

def _get_target_mode(args):
    """Construct TargetMode based on CLI arguments."""
    target_str = getattr(args, 'target', None) or 'kind16'
    return load_target(target_str)

def _parser_args(args):
    """Extract parser/parser_cmd from CLI args."""
    parser = getattr(args, 'parser', None)
    parser_cmd = getattr(args, 'parser_cmd', None)
    return parser, parser_cmd


def cmd_migrate(args):
    """Run the migration step."""
    parser, parser_cmd = _parser_args(args)
    target = _get_target_mode(args)
    run_migration(
        recipe_path=args.recipe,
        output_dir=args.output_dir,
        target_mode=target,
        dry_run=args.dry_run,
        project_root=args.project_root,
        parser=parser,
        parser_cmd=parser_cmd,
    )


def cmd_diverge(args):
    """Report every co-family pair whose migrated text differs."""
    parser, parser_cmd = _parser_args(args)
    target = _get_target_mode(args)
    report = run_divergence_report(
        recipe_path=args.recipe,
        target_mode=target,
        project_root=args.project_root,
        parser=parser,
        parser_cmd=parser_cmd,
    )
    total = len(report)
    # Optional filtering on diff content.
    if args.grep:
        pat = re.compile(args.grep, re.IGNORECASE)
        report = [r for r in report if any(pat.search(l) for l in r['diff'])]
    if args.exclude:
        pat = re.compile(args.exclude, re.IGNORECASE)
        report = [r for r in report if not any(pat.search(l) for l in r['diff'])]

    for entry in report:
        header = (f'### {entry["other"]} vs {entry["canonical"]}'
                  f' → {entry["target"]} (+{len(entry["diff"])})')
        print(header)
        diff = entry['diff'] if args.full else entry['diff'][:args.context]
        for line in diff:
            print(line[:args.max_width])
        if not args.full and len(entry['diff']) > args.context:
            print(f'  ...{len(entry["diff"]) - args.context} more')
        print()

    shown = len(report)
    if args.grep or args.exclude:
        print(f'{shown} shown / {total} divergent pairs')
    else:
        print(f'{total} divergent pairs')


def cmd_converge(args):
    """Verify that migrated files on disk converge."""
    parser, parser_cmd = _parser_args(args)
    target = _get_target_mode(args)
    report = run_convergence_report(
        recipe_path=args.recipe,
        output_dir=args.output_dir,
        target_mode=target,
        project_root=args.project_root,
        parser=parser,
        parser_cmd=parser_cmd,
    )
    total = len(report)
    if args.grep:
        pat = re.compile(args.grep, re.IGNORECASE)
        report = [r for r in report if any(pat.search(l) for l in r['diff'])]
    if args.exclude:
        pat = re.compile(args.exclude, re.IGNORECASE)
        report = [r for r in report if not any(pat.search(l) for l in r['diff'])]

    for entry in report:
        if entry['status'] == 'missing':
            print(f'### MISSING {entry["target"]} '
                  f'(expected from {entry["canonical"]})')
            print()
            continue
        header = (f'### {entry["other"]} vs {entry["canonical"]}'
                  f' → {entry["target"]} (+{len(entry["diff"])})')
        print(header)
        diff = entry['diff'] if args.full else entry['diff'][:args.context]
        for line in diff:
            print(line[:args.max_width])
        if not args.full and len(entry['diff']) > args.context:
            print(f'  ...{len(entry["diff"]) - args.context} more')
        print()

    diverged = sum(1 for r in report if r['status'] == 'diverged')
    missing = sum(1 for r in report if r['status'] == 'missing')
    shown = len(report)
    summary_total = sum(1 for _ in range(total))  # == total
    if args.grep or args.exclude:
        print(f'{shown} shown / {summary_total} entries '
              f'({diverged} diverged, {missing} missing on disk)')
    else:
        print(f'{diverged} diverged, {missing} missing on disk')


def _is_fixed_form_comment(line: str) -> bool:
    """A fixed-form line is a comment if its first character is C/c/*/!
    OR if its first non-whitespace character is ``!`` (the inline-comment
    marker can also start a whole-line comment when it stands alone)."""
    if not line:
        return True
    if line[0] in ('C', 'c', '*', '!'):
        return True
    stripped = line.lstrip()
    return stripped.startswith('!')


def _is_free_form_comment(line: str) -> bool:
    return not line or line.lstrip().startswith('!')


def cmd_verify(args):
    """Verify migrated output: residual types, literals, column width."""
    out_dir = args.output_dir
    src_dir = out_dir / 'src'
    if not src_dir.is_dir():
        # Fall back to flat layout
        src_dir = out_dir
    errors = 0

    target_mode = _get_target_mode(args)

    f_files = sorted(src_dir.glob('*.f'))
    f90_files = sorted(list(src_dir.glob('*.f90')) + list(src_dir.glob('*.F90')))
    all_files = f_files + f90_files

    print(f'Verifying {len(all_files)} files in {src_dir}')
    print()

    # Determine if the target uses module-based constructors
    is_constructor_based = target_mode.real_constructor is not None

    # Check residual precision types in code lines
    print('Residual precision types (code lines):')
    residuals = 0
    for f in f_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if _is_fixed_form_comment(line):
                continue
            if re.search(r'DOUBLE\s+PRECISION|COMPLEX\*16|COMPLEX\*8|DOUBLE\s+COMPLEX|REAL\*[48]',
                         line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                residuals += 1
    for f in f90_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if _is_free_form_comment(line):
                continue
            if re.search(r'kind\s*\(\s*1\.[de]0\s*\)', line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                residuals += 1
            if re.search(r'\bdouble\s+precision\b', line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                residuals += 1
    if residuals == 0:
        print('  OK')
    else:
        errors += residuals
    print()

    # Build constructor-stripping patterns from target_mode
    def _strip_constructors(line: str) -> str:
        """Remove target-type constructor wrappers from a line."""
        if not target_mode.real_constructor:
            return line
        ctor = re.escape(target_mode.real_constructor)
        line = re.sub(rf"{ctor}\(\s*'[^']*'\s*\)", '', line, flags=re.IGNORECASE)
        line = re.sub(rf'{ctor}\(\s*[^)]*\s*\)', '', line, flags=re.IGNORECASE)
        if target_mode.complex_constructor:
            cctor = re.escape(target_mode.complex_constructor)
            line = re.sub(rf"{cctor}\(\s*'[^']*'\s*\)", '', line, flags=re.IGNORECASE)
            line = re.sub(rf'{cctor}\(\s*[^)]*\s*\)', '', line, flags=re.IGNORECASE)
        return line

    # Check residual D-exponent literals in code lines
    print('Residual D-exponent literals (code lines):')
    d_lits = 0
    for f in f_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if _is_fixed_form_comment(line):
                continue
            cleaned_line = _strip_constructors(line)
            if re.search(r'[0-9]\.[0-9]*[Dd][+-]?[0-9]', cleaned_line):
                print(f'  {f.name}:{i}: {line.strip()}')
                d_lits += 1
    for f in f90_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if _is_free_form_comment(line):
                continue
            cleaned_line = _strip_constructors(line)
            # In constructor mode, also reject _wp suffixed literals
            if is_constructor_based and re.search(r'\d+\.\d*_wp', cleaned_line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                d_lits += 1
    if d_lits == 0:
        print('  OK')
    else:
        errors += d_lits
    print()

    if is_constructor_based:
        # Constructor-based target: residual unconverted FP PARAMETER /
        # DATA statements (those that mention a value-shaped numeric and
        # are not commented out).
        print('Residual FP PARAMETER/DATA (code lines):')
        leftover = 0
        for f in all_files:
            text = f.read_text(errors='replace')
            is_fixed = f.suffix.lower() == '.f'
            for i, line in enumerate(text.splitlines(), 1):
                if is_fixed and _is_fixed_form_comment(line):
                    continue
                if not is_fixed and _is_free_form_comment(line):
                    continue
                m = re.match(r'\s+(PARAMETER|DATA)\b', line, re.IGNORECASE)
                if not m:
                    continue
                cleaned = _strip_constructors(line)
                if re.search(r'\d+\.\d*[DdEe][+-]?\d+|\d*\.\d+', cleaned):
                    print(f'  {f.name}:{i}: {line.strip()}')
                    leftover += 1
        if leftover == 0:
            print('  OK')
        else:
            errors += leftover
        print()

    # Check column width for fixed-form code lines
    print('Column overflow (code lines > 72 chars):')
    overflows = 0
    for f in f_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if _is_fixed_form_comment(line):
                continue
            if len(line) > 72:
                print(f'  {f.name}:{i}: {len(line)} chars')
                overflows += 1
    if overflows == 0:
        print('  OK')
    else:
        errors += overflows
    print()

    if errors:
        print(f'FAILED: {errors} issue(s)')
        sys.exit(1)
    else:
        print('PASSED')


def _generate_cmake(output_dir: Path, lib_name: str, target_mode,
                    common_files: list[str], precision_files: list[str],
                    language: str = 'fortran',
                    project_root: Path | None = None):
    """Generate a self-contained CMakeLists.txt in the output directory."""
    pmap = target_mode.prefix_map
    real_pfx = pmap['R'].lower()
    precision_lib = f'{real_pfx}{lib_name}'
    common_lib = f'{lib_name}_common'

    common_list = '\n    '.join(sorted(common_files))
    precision_list = '\n    '.join(sorted(precision_files))

    if language == 'c':
        cmake = f"""\
cmake_minimum_required(VERSION 3.20)
project({precision_lib} C)

# --- Compiler flags ---
set(CMAKE_C_FLAGS "${{CMAKE_C_FLAGS}} -w")

# --- MPI (optional, for libraries like BLACS) ---
find_package(MPI COMPONENTS C QUIET)
if(MPI_C_FOUND)
    include_directories(${{MPI_C_INCLUDE_DIRS}})
endif()

# --- Common (type-independent) library ---
set(COMMON_SOURCES
    {common_list}
)

# --- Precision-specific library ---
set(PRECISION_SOURCES
    {precision_list}
)

# Header include path
include_directories(${{CMAKE_CURRENT_SOURCE_DIR}}/src)

if(COMMON_SOURCES)
    add_library({common_lib} STATIC ${{COMMON_SOURCES}})
endif()

add_library({precision_lib} STATIC ${{PRECISION_SOURCES}})
if(TARGET {common_lib})
    target_link_libraries({precision_lib} PUBLIC {common_lib})
endif()

# --- Install rules ---
install(TARGETS {precision_lib} ARCHIVE DESTINATION lib)
if(TARGET {common_lib})
    install(TARGETS {common_lib} ARCHIVE DESTINATION lib)
endif()
"""
    else:
        # If multifloats, we need to link against the multifloats library
        # AND build the la_constants_mf / la_xisnan_mf helper modules
        # that the migrated source depends on for la_constants USE clauses.
        mf_link = ""
        mf_deps = ""
        if target_mode.module_name is not None:
            # Resolve absolute paths to external dependencies so the
            # generated CMakeLists.txt works from any output directory.
            _root = project_root or Path.cwd()
            _mf_default = str((_root / 'external' / 'multifloats').resolve())
            _helpers_default = str((_root / 'external' / 'lapack-3.12.1' / 'SRC').resolve())
            mf_link = f"""
# Find or link multifloats. The location can be overridden via the
# MULTIFLOATS_DIR cache variable; otherwise we use the resolved path
# from the project root at generation time.
if(NOT DEFINED MULTIFLOATS_DIR)
    set(MULTIFLOATS_DIR "{_mf_default}"
        CACHE PATH "Path to the external multifloats source tree")
endif()
if(EXISTS "${{MULTIFLOATS_DIR}}/CMakeLists.txt")
    add_subdirectory(${{MULTIFLOATS_DIR}} multifloats_build)
else()
    message(FATAL_ERROR
        "multifloats library not found at ${{MULTIFLOATS_DIR}}. Set "
        "-DMULTIFLOATS_DIR=/path/to/multifloats to override.")
endif()

# Build the la_constants_mf and la_xisnan_mf helper modules. These
# re-export multifloats's MF_* constants under the DD/ZZ-prefixed names
# that the migrated LAPACK source uses via its rewritten
# ``USE LA_CONSTANTS_MF`` clause.
set(MF_HELPERS_DIR "{_helpers_default}"
    CACHE PATH "Directory containing la_constants_mf.f90 / la_xisnan_mf.f90")
if(EXISTS "${{MF_HELPERS_DIR}}/la_constants_mf.f90")
    add_library(la_constants_mf STATIC
        "${{MF_HELPERS_DIR}}/la_constants_mf.f90")
    set_target_properties(la_constants_mf PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories(la_constants_mf PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    if(TARGET multifloats)
        target_link_libraries(la_constants_mf PUBLIC multifloats)
    endif()
endif()
if(EXISTS "${{MF_HELPERS_DIR}}/la_xisnan_mf.f90")
    add_library(la_xisnan_mf STATIC
        "${{MF_HELPERS_DIR}}/la_xisnan_mf.f90")
    set_target_properties(la_xisnan_mf PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories(la_xisnan_mf PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    if(TARGET multifloats)
        target_link_libraries(la_xisnan_mf PUBLIC multifloats)
    endif()
endif()
"""
            mf_deps = f"""
if(TARGET multifloats)
    target_link_libraries({precision_lib} PUBLIC multifloats)
endif()
if(TARGET la_constants_mf)
    target_link_libraries({precision_lib} PUBLIC la_constants_mf)
endif()
if(TARGET la_xisnan_mf)
    target_link_libraries({precision_lib} PUBLIC la_xisnan_mf)
endif()
"""

        cmake = f"""\
cmake_minimum_required(VERSION 3.20)
project({precision_lib} Fortran)

# --- Compiler flags ---
set(CMAKE_Fortran_FLAGS "${{CMAKE_Fortran_FLAGS}} -w")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${{CMAKE_Fortran_FLAGS}} -std=legacy")
endif()

# Enable Fortran preprocessing for .F90 files
set(CMAKE_Fortran_PREPROCESS ON)

# --- MPI (optional, for libraries like MUMPS that INCLUDE 'mpif.h') ---
find_package(MPI COMPONENTS Fortran QUIET)

# Detect 80-bit extended precision (KIND=10) support
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    include(CheckFortranSourceCompiles)
    check_fortran_source_compiles("
        program test_real10
        integer, parameter :: ep = selected_real_kind(18, 4931)
        real(ep) :: x
        x = 1.0_ep
        end program
    " HAVE_REAL10 SRC_EXT F90)
    if(HAVE_REAL10)
        add_compile_definitions(HAVE_REAL10)
    endif()
endif()
{mf_link}
# --- Common (type-independent) library ---
set(COMMON_SOURCES
    {common_list}
)

# --- Precision-specific library ---
set(PRECISION_SOURCES
    {precision_list}
)

if(COMMON_SOURCES)
    add_library({common_lib} STATIC ${{COMMON_SOURCES}})
    set_target_properties({common_lib} PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories({common_lib} PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
endif()

add_library({precision_lib} STATIC ${{PRECISION_SOURCES}})
set_target_properties({precision_lib} PROPERTIES
    Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
target_include_directories({precision_lib} PUBLIC
    $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>
    $<BUILD_INTERFACE:${{CMAKE_CURRENT_SOURCE_DIR}}>)
if(TARGET {common_lib})
    target_link_libraries({precision_lib} PUBLIC {common_lib})
endif()
if(MPI_Fortran_FOUND)
    target_link_libraries({precision_lib} PUBLIC MPI::MPI_Fortran)
    if(TARGET {common_lib})
        target_link_libraries({common_lib} PUBLIC MPI::MPI_Fortran)
    endif()
endif()
{mf_deps}

# --- Install rules ---
install(TARGETS {precision_lib} ARCHIVE DESTINATION lib)
if(TARGET {common_lib})
    install(TARGETS {common_lib} ARCHIVE DESTINATION lib)
endif()
"""
    (output_dir / 'CMakeLists.txt').write_text(cmake)


def cmd_build(args):
    """Generate CMake project and build static libraries."""
    output_dir = args.output_dir
    target_mode = _get_target_mode(args)
    src_dir = output_dir / 'src'
    if not src_dir.is_dir():
        src_dir = output_dir

    config = load_recipe(args.recipe, args.project_root)
    lib_name = config.library

    # Classify source files into common vs precision-specific
    symbols = scan_symbols(config.source_dir, config.language,
                           config.extensions, config.library_path,
                           extra_c_return_types=tuple(config.c_return_types))
    classification = classify_symbols(symbols)
    independent = classification.independent

    if config.language == 'c':
        files = sorted(list(src_dir.glob('*.c')))
    else:
        # Honor the recipe's extensions list (normalized to lowercase in
        # load_recipe) so libraries that use .F (MUMPS) or any
        # non-default extension are picked up. Case-insensitive match on
        # the actual filename suffix.
        allowed = {e.lower() for e in config.extensions}
        files = sorted(
            p for p in src_dir.iterdir()
            if p.is_file() and p.suffix.lower() in allowed
        )
    common_files, precision_files = [], []
    for f in files:
        rel = f.relative_to(output_dir)
        if f.stem.upper() in independent:
            common_files.append(str(rel))
        else:
            precision_files.append(str(rel))

    print(f'Generating CMake project in {output_dir}/')
    print(f'  Common:    {len(common_files)} files')
    print(f'  Precision: {len(precision_files)} files')

    proj_root = (args.project_root or args.recipe.resolve().parent.parent)
    _generate_cmake(output_dir, lib_name, target_mode, common_files, precision_files,
                    language=config.language, project_root=proj_root)

    # Configure and build
    build_dir = output_dir / '_build'
    cmake_cmd = 'cmake'

    # Configure
    configure_args = [
        cmake_cmd, '-S', str(output_dir), '-B', str(build_dir),
        '-DCMAKE_BUILD_TYPE=Release',
    ]
    if config.language != 'c' and args.fc:
        configure_args.append(f'-DCMAKE_Fortran_COMPILER={args.fc}')

    print(f'\nConfiguring...')
    r = subprocess.run(configure_args, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        sys.exit(1)

    # Build (parallel)
    jobs = os.cpu_count() or 4
    print(f'Building ({jobs} parallel jobs)...')
    r = subprocess.run(
        [cmake_cmd, '--build', str(build_dir), '-j', str(jobs)],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        # Show failures
        for line in r.stderr.splitlines():
            if 'error' in line.lower() or 'Error' in line:
                print(f'  {line}')
        # Count pass/fail from build output
        lines = r.stdout.splitlines() + r.stderr.splitlines()
        errors = [l for l in lines if 'Error:' in l]
        print(f'\nBuild failed with {len(errors)} error(s)')
        print(f'Full log: {build_dir}/')
        sys.exit(1)

    # Report results
    pmap = target_mode.prefix_map
    real_pfx = pmap['R'].lower()
    precision_lib_name = f'lib{real_pfx}{lib_name}.a'
    common_lib_name = f'lib{lib_name}_common.a'

    print(f'\nBuild succeeded:')
    for name in [common_lib_name, precision_lib_name]:
        matches = list(build_dir.rglob(name))
        if matches:
            p = matches[0]
            size = p.stat().st_size
            print(f'  {name}: {size // 1024}K')

    print(f'\nLibraries in {build_dir}/')


def cmd_run(args):
    """Run the full pipeline: migrate → converge → verify → build."""
    work_dir = args.work_dir
    output_dir = work_dir / 'output'
    src_dir = output_dir / 'src'
    src_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 60)
    print('  Step 1: Migrate')
    print('=' * 60)
    args.output_dir = src_dir
    args.dry_run = False
    cmd_migrate(args)

    print()
    print('=' * 60)
    print('  Step 2: Convergence')
    print('=' * 60)
    args.output_dir = src_dir
    # Set defaults for convergence report display options
    if not hasattr(args, 'grep'):
        args.grep = None
    if not hasattr(args, 'exclude'):
        args.exclude = None
    if not hasattr(args, 'context'):
        args.context = 8
    if not hasattr(args, 'full'):
        args.full = False
    if not hasattr(args, 'max_width'):
        args.max_width = 200
    cmd_converge(args)

    print()
    print('=' * 60)
    print('  Step 3: Verify')
    print('=' * 60)
    args.output_dir = output_dir
    try:
        cmd_verify(args)
    except SystemExit as e:
        if e.code:
            print('Verify failed, continuing...')

    print()
    print('=' * 60)
    print('  Step 4: Build (CMake)')
    print('=' * 60)
    args.output_dir = output_dir
    cmd_build(args)


# Topologically sorted library build order for the unified CMake project.
# Each entry is (library_name, recipe_filename).
LIBRARY_ORDER = [
    ('blas',        'blas.yaml'),
    ('blacs',       'blacs.yaml'),
    ('lapack',      'lapack.yaml'),
    ('ptzblas',     'ptzblas.yaml'),
    ('pbblas',      'pbblas.yaml'),
    ('pblas',       'pblas.yaml'),
    ('scalapack',   'scalapack.yaml'),
    ('scalapack_c', 'scalapack_c.yaml'),
]


def cmd_stage(args):
    """Migrate all libraries into a structured staging directory.

    Produces a self-contained directory that can be built with:
        cmake -S <staging> -B <staging>/build && cmake --build <staging>/build -j
    """
    staging_dir = args.staging_dir.resolve()
    target_mode = _get_target_mode(args)
    parser, parser_cmd = _parser_args(args)
    proj_root = (args.project_root or Path(__file__).resolve().parent.parent.parent)
    recipes_dir = proj_root / 'recipes'

    # Determine which libraries to stage
    if args.libraries:
        lib_set = set(args.libraries)
        libraries = [(n, r) for n, r in LIBRARY_ORDER if n in lib_set]
    else:
        libraries = list(LIBRARY_ORDER)

    staged = []
    for lib_name, recipe_file in libraries:
        recipe_path = recipes_dir / recipe_file
        if not recipe_path.exists():
            print(f'Warning: recipe {recipe_path} not found, skipping {lib_name}')
            continue

        lib_dir = staging_dir / lib_name
        src_dir = lib_dir / 'src'
        src_dir.mkdir(parents=True, exist_ok=True)

        print(f'\n{"=" * 60}')
        print(f'  Migrating: {lib_name}')
        print(f'{"=" * 60}')

        # Run migration
        run_migration(
            recipe_path=recipe_path,
            output_dir=src_dir,
            target_mode=target_mode,
            dry_run=False,
            project_root=proj_root,
            parser=parser,
            parser_cmd=parser_cmd,
        )

        # Classify files into common vs precision-specific
        config = load_recipe(recipe_path, proj_root)
        symbols = scan_symbols(config.source_dir, config.language,
                               config.extensions, config.library_path,
                               extra_c_return_types=tuple(config.c_return_types))
        classification = classify_symbols(symbols)
        independent = classification.independent

        if config.language == 'c':
            files = sorted(list(src_dir.glob('*.c')))
        else:
            files = sorted(
                list(src_dir.glob('*.f')) + list(src_dir.glob('*.f90'))
                + list(src_dir.glob('*.F90'))
            )

        # MF helper modules are built as separate targets; exclude from manifests.
        _mf_helpers = {'la_constants_mf', 'la_xisnan_mf'}
        # la_constants_ep provides extended/quad precision constants —
        # not needed by multifloats target, and may fail on compilers
        # that don't support REAL(KIND=16).
        if target_mode.module_name is not None:
            _mf_helpers.update({'la_constants_ep', 'la_xisnan_ep'})

        common_files, precision_files = [], []
        for f in files:
            if f.stem in _mf_helpers:
                continue
            rel = f'src/{f.name}'
            if f.stem.upper() in independent:
                common_files.append(rel)
            else:
                precision_files.append(rel)

        # Write manifest.cmake
        common_list = '\n    '.join(common_files) if common_files else ''
        precision_list = '\n    '.join(precision_files) if precision_files else ''
        manifest = f"""\
set({lib_name}_COMMON_SOURCES
    {common_list}
)

set({lib_name}_PRECISION_SOURCES
    {precision_list}
)

set({lib_name}_LANGUAGE {config.language})
"""
        (lib_dir / 'manifest.cmake').write_text(manifest)
        print(f'  Manifest: {len(common_files)} common, '
              f'{len(precision_files)} precision files')
        staged.append(lib_name)

    # Copy MF helper modules into staging so it's self-contained
    pmap = target_mode.prefix_map
    lib_prefix = pmap['R'].lower()
    needs_mf = target_mode.module_name is not None
    staged_list = ';'.join(staged)

    helpers_src = proj_root / 'external' / 'lapack-3.12.1' / 'SRC'
    helpers_dst = staging_dir / '_helpers'
    if needs_mf:
        helpers_dst.mkdir(exist_ok=True)
        for name in ['la_constants_mf.f90', 'la_xisnan_mf.f90']:
            src = helpers_src / name
            if src.exists():
                shutil.copy2(src, helpers_dst / name)
        # Copy multifloats bridge files (C++ bridge header + MPI registration)
        mf_local = proj_root / 'external' / 'multifloats'
        bridge_h = mf_local / 'include' / 'multifloats_bridge.h'
        mpi_cpp = mf_local / 'src' / 'multifloats_mpi.cpp'
        if bridge_h.exists():
            shutil.copy2(bridge_h, helpers_dst / bridge_h.name)
        if mpi_cpp.exists():
            shutil.copy2(mpi_cpp, helpers_dst / mpi_cpp.name)

    target_config = f"""\
# Generated by: python -m pyengine stage --target {target_mode.name}
set(TARGET_NAME "{target_mode.name}")
set(LIB_PREFIX "{lib_prefix}")
set(NEEDS_MULTIFLOATS {'TRUE' if needs_mf else 'FALSE'})
set(C_AS_CXX {'TRUE' if needs_mf else 'FALSE'})
set(MF_HELPERS_DIR "${{CMAKE_CURRENT_SOURCE_DIR}}/_helpers")
set(STAGED_LIBRARIES {staged_list})
"""
    (staging_dir / 'target_config.cmake').write_text(target_config)

    # Copy CMake files to staging directory
    cmake_dir = proj_root / 'cmake'
    for cmake_file in ['CMakeLists.txt', 'FortranCompiler.cmake']:
        src = cmake_dir / cmake_file
        if src.exists():
            shutil.copy2(src, staging_dir / cmake_file)
        else:
            print(f'Warning: {src} not found')

    print(f'\n{"=" * 60}')
    print(f'  Staging complete: {len(staged)} libraries')
    print(f'{"=" * 60}')
    print(f'  Target:  {target_mode.name} (prefix: {lib_prefix})')
    print(f'  Output:  {staging_dir}')
    print(f'\nTo build:')
    print(f'  cmake -S {staging_dir} -B {staging_dir}/build -DCMAKE_BUILD_TYPE=Release')
    print(f'  cmake --build {staging_dir}/build -j')


def _add_parser_args(p):
    """Add --parser and --parser-cmd arguments to a subparser."""
    p.add_argument(
        '--parser', default=None,
        choices=['flang', 'gfortran'],
        help='Parse tree backend for Fortran migration '
             '(default: regex-only, no compiler)')
    p.add_argument(
        '--parser-cmd', default=None,
        help='Explicit path to the parser compiler binary '
             '(overrides PATH lookup)')

def _add_target_args(p):
    p.add_argument('--target', type=str, default='kind16',
                   help='Target name (e.g. "multifloats", "kind16") or path to a target .yaml file')

def main():
    parser = argparse.ArgumentParser(
        prog='pyengine',
        description='General-purpose type migration for numerical libraries'
    )
    sub = parser.add_subparsers(dest='command', required=True)

    # --- migrate ---
    p = sub.add_parser('migrate', help='Migrate source files')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('output_dir', type=Path, help='Output directory')
    _add_target_args(p)
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--project-root', type=Path, default=None)
    _add_parser_args(p)
    p.set_defaults(func=cmd_migrate)

    # --- diverge ---
    p = sub.add_parser('diverge',
                       help='Report co-family pairs with differing output')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    _add_target_args(p)
    p.add_argument('--project-root', type=Path, default=None)
    p.add_argument('--grep', default=None,
                   help='Regex: only show entries with diff matching')
    p.add_argument('--exclude', default=None,
                   help='Regex: drop entries whose diff matches')
    p.add_argument('--context', type=int, default=8,
                   help='Max diff lines per entry (default 8)')
    p.add_argument('--full', action='store_true',
                   help='Print full diff per entry (ignores --context)')
    p.add_argument('--max-width', type=int, default=200,
                   help='Truncate each diff line to this many chars')
    _add_parser_args(p)
    p.set_defaults(func=cmd_diverge)

    # --- converge ---
    p = sub.add_parser('converge',
                       help='Post-migration verification against on-disk '
                            'files with a light whitespace-only normalizer')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('output_dir', type=Path,
                   help='Directory holding migrated output')
    _add_target_args(p)
    p.add_argument('--project-root', type=Path, default=None)
    p.add_argument('--grep', default=None,
                   help='Regex: only show entries with diff matching')
    p.add_argument('--exclude', default=None,
                   help='Regex: drop entries whose diff matches')
    p.add_argument('--context', type=int, default=8,
                   help='Max diff lines per entry (default 8)')
    p.add_argument('--full', action='store_true',
                   help='Print full diff per entry (ignores --context)')
    p.add_argument('--max-width', type=int, default=200,
                   help='Truncate each diff line to this many chars')
    _add_parser_args(p)
    p.set_defaults(func=cmd_converge)

    # --- verify ---
    p = sub.add_parser('verify', help='Verify migrated output')
    p.add_argument('output_dir', type=Path)
    p.set_defaults(func=cmd_verify)

    # --- build ---
    p = sub.add_parser('build', help='Generate CMake project and build')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('output_dir', type=Path)
    _add_target_args(p)
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--project-root', type=Path, default=None)
    p.set_defaults(func=cmd_build)

    # --- run (full pipeline) ---
    p = sub.add_parser('run', help='Run full pipeline')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('work_dir', type=Path, help='Working directory')
    _add_target_args(p)
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--project-root', type=Path, default=None)
    _add_parser_args(p)
    p.set_defaults(func=cmd_run)

    # --- stage (migrate all + unified cmake) ---
    p = sub.add_parser('stage',
                       help='Migrate all libraries into a staging directory '
                            'with a unified CMake build')
    p.add_argument('staging_dir', type=Path,
                   help='Output staging directory')
    _add_target_args(p)
    p.add_argument('--project-root', type=Path, default=None)
    p.add_argument('--libraries', nargs='+', default=None,
                   help='Subset of libraries to migrate (default: all)')
    _add_parser_args(p)
    p.set_defaults(func=cmd_stage)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
