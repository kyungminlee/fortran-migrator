"""CLI entry point for the general-purpose migration pipeline.

Usage:
    uv run python -m pyengine migrate  recipes/blas.yaml output/ --kind 16
    uv run python -m pyengine verify   output/
    uv run python -m pyengine build    recipes/blas.yaml output/ --kind 16
    uv run python -m pyengine run      recipes/blas.yaml work/ --kind 16
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

from .config import load_recipe
from .pipeline import (
    run_convergence_report, run_divergence_report, run_migration,
)
from .prefix_classifier import classify_symbols, PREFIX_MAP
from .symbol_scanner import scan_symbols


def _parser_args(args):
    """Extract parser/parser_cmd from CLI args."""
    parser = getattr(args, 'parser', None)
    parser_cmd = getattr(args, 'parser_cmd', None)
    return parser, parser_cmd


def cmd_migrate(args):
    """Run the migration step."""
    parser, parser_cmd = _parser_args(args)
    run_migration(
        recipe_path=args.recipe,
        output_dir=args.output_dir,
        target_kind=args.kind,
        dry_run=args.dry_run,
        project_root=args.project_root,
        parser=parser,
        parser_cmd=parser_cmd,
    )


def cmd_diverge(args):
    """Report every co-family pair whose migrated text differs.

    Shows one entry per divergent pair with the +/- lines from the
    canonicalized unified diff (canonical = D/Z half, kept on disk;
    other = S/C half, checked against it).
    """
    parser, parser_cmd = _parser_args(args)
    report = run_divergence_report(
        recipe_path=args.recipe,
        target_kind=args.kind,
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
    """Verify that migrated files on disk converge with a fresh
    migration of their S/C co-family siblings under a light
    normalizer (whitespace/case only).

    This is the authoritative post-migration check: the D/Z
    canonical is read from disk, the S/C sibling is re-migrated
    in memory, and the two are compared without the heavy prefix
    collapse / literal rewrites that ``diverge`` uses to mask
    pre-existing BLAS/LAPACK source drift. Any entry here means
    either the migrator left precision-specific material behind or
    a genuine algorithmic difference exists between the halves.
    """
    parser, parser_cmd = _parser_args(args)
    report = run_convergence_report(
        recipe_path=args.recipe,
        output_dir=args.output_dir,
        target_kind=args.kind,
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


def cmd_verify(args):
    """Verify migrated output: residual types, literals, column width."""
    out_dir = args.output_dir
    src_dir = out_dir / 'src'
    if not src_dir.is_dir():
        # Fall back to flat layout
        src_dir = out_dir
    errors = 0

    f_files = sorted(src_dir.glob('*.f'))
    f90_files = sorted(list(src_dir.glob('*.f90')) + list(src_dir.glob('*.F90')))
    all_files = f_files + f90_files

    print(f'Verifying {len(all_files)} files in {src_dir}')
    print()

    # Check residual precision types in code lines
    print('Residual precision types (code lines):')
    residuals = 0
    for f in f_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if not line or line[0] in ('C', 'c', '*', '!'):
                continue
            if re.search(r'DOUBLE\s+PRECISION|COMPLEX\*16|COMPLEX\*8|DOUBLE\s+COMPLEX',
                         line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                residuals += 1
    for f in f90_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if re.search(r'kind\s*\(\s*1\.[de]0\s*\)', line, re.IGNORECASE):
                print(f'  {f.name}:{i}: {line.strip()}')
                residuals += 1
    if residuals == 0:
        print('  OK')
    else:
        errors += residuals
    print()

    # Check residual D-exponent literals in code lines
    print('Residual D-exponent literals (code lines):')
    d_lits = 0
    for f in f_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if not line or line[0] in ('C', 'c', '*', '!'):
                continue
            if re.search(r'[0-9]\.[0-9]*[Dd][+-]?[0-9]', line):
                print(f'  {f.name}:{i}: {line.strip()}')
                d_lits += 1
    if d_lits == 0:
        print('  OK')
    else:
        errors += d_lits
    print()

    # Check column width for fixed-form code lines
    print('Column overflow (code lines > 72 chars):')
    overflows = 0
    for f in f_files:
        for i, line in enumerate(f.read_text(errors='replace').splitlines(), 1):
            if not line or line[0] in ('C', 'c', '*', '!'):
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


def _generate_cmake(output_dir: Path, lib_name: str, kind: int,
                    common_files: list[str], precision_files: list[str]):
    """Generate a self-contained CMakeLists.txt in the output directory."""
    pmap = PREFIX_MAP[kind]
    real_pfx = pmap['R'].lower()
    precision_lib = f'{real_pfx}{lib_name}'
    common_lib = f'{lib_name}_common'

    common_list = '\n    '.join(sorted(common_files))
    precision_list = '\n    '.join(sorted(precision_files))

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
    $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
if(TARGET {common_lib})
    target_link_libraries({precision_lib} PUBLIC {common_lib})
endif()

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
    kind = args.kind
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

    files = sorted(
        list(src_dir.glob('*.f')) + list(src_dir.glob('*.f90'))
        + list(src_dir.glob('*.F90'))
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

    _generate_cmake(output_dir, lib_name, kind, common_files, precision_files)

    # Configure and build
    build_dir = output_dir / '_build'
    cmake_cmd = 'cmake'

    # Configure
    configure_args = [
        cmake_cmd, '-S', str(output_dir), '-B', str(build_dir),
        '-DCMAKE_BUILD_TYPE=Release',
    ]
    if args.fc:
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
    pmap = PREFIX_MAP[kind]
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
    """Run the full pipeline: migrate → verify → build."""
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
    print('  Step 2: Verify')
    print('=' * 60)
    args.output_dir = output_dir
    try:
        cmd_verify(args)
    except SystemExit as e:
        if e.code:
            print('Verify failed, continuing...')

    print()
    print('=' * 60)
    print('  Step 3: Build (CMake)')
    print('=' * 60)
    args.output_dir = output_dir
    cmd_build(args)


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
    p.add_argument('--kind', type=int, default=16, choices=[10, 16])
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--project-root', type=Path, default=None)
    _add_parser_args(p)
    p.set_defaults(func=cmd_migrate)

    # --- diverge ---
    p = sub.add_parser('diverge',
                       help='Report co-family pairs with differing output')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('--kind', type=int, default=16, choices=[10, 16])
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
    p.add_argument('--kind', type=int, default=16, choices=[10, 16])
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
    p.add_argument('--kind', type=int, default=16, choices=[10, 16])
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--project-root', type=Path, default=None)
    p.set_defaults(func=cmd_build)

    # --- run (full pipeline) ---
    p = sub.add_parser('run', help='Run full pipeline')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('work_dir', type=Path, help='Working directory')
    p.add_argument('--kind', type=int, default=16, choices=[10, 16])
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--project-root', type=Path, default=None)
    _add_parser_args(p)
    p.set_defaults(func=cmd_run)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
