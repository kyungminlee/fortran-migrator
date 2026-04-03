"""CLI entry point for the general-purpose migration pipeline.

Usage:
    uv run python -m pyengine migrate  recipes/blas.yaml output/ --kind 16
    uv run python -m pyengine verify   output/
    uv run python -m pyengine compile  output/ [--fc gfortran]
    uv run python -m pyengine build    output/ --lib-name libqblas.a
    uv run python -m pyengine run      recipes/blas.yaml work/ --kind 16
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from .config import load_recipe
from .pipeline import run_migration
from .prefix_classifier import classify_symbols, build_rename_map, PREFIX_MAP
from .symbol_scanner import scan_symbols


def cmd_migrate(args):
    """Run the migration step."""
    run_migration(
        recipe_path=args.recipe,
        output_dir=args.output_dir,
        target_kind=args.kind,
        dry_run=args.dry_run,
        project_root=args.project_root,
    )


def cmd_verify(args):
    """Verify migrated output: residual types, literals, column width."""
    out_dir = args.output_dir
    errors = 0

    if not out_dir.is_dir():
        print(f'Error: output directory not found: {out_dir}', file=sys.stderr)
        sys.exit(1)

    f_files = sorted(out_dir.glob('*.f'))
    f90_files = sorted(out_dir.glob('*.f90'))
    all_files = f_files + f90_files

    print(f'Verifying {len(all_files)} files in {out_dir}')
    print()

    # Check residual precision types in code lines
    print('Residual precision types (code lines):')
    residuals = 0
    import re
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


def cmd_compile(args):
    """Compile all migrated files to verify syntax."""
    out_dir = args.output_dir
    fc = args.fc
    fflags = ['-c', '-O0', '-w']

    if not shutil.which(fc):
        print(f'Error: Fortran compiler not found: {fc}', file=sys.stderr)
        sys.exit(1)

    files = sorted(list(out_dir.glob('*.f')) + list(out_dir.glob('*.f90')))
    print(f'Compiling {len(files)} files with {fc}')

    build_dir = out_dir / '_build'
    build_dir.mkdir(exist_ok=True)

    passed, failed = 0, 0
    for src in files:
        obj = build_dir / (src.stem + '.o')
        result = subprocess.run(
            [fc] + fflags + ['-o', str(obj), str(src)],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            passed += 1
        else:
            failed += 1
            print(f'  FAIL: {src.name}')
            for line in result.stderr.splitlines()[:5]:
                print(f'    {line}')

    print(f'\nPASS: {passed}  FAIL: {failed}')
    if failed:
        sys.exit(1)


def cmd_build(args):
    """Build static libraries (common + precision-specific)."""
    out_dir = args.output_dir
    fc = args.fc
    ar = args.ar
    fflags = ['-c', '-O2', '-w']
    kind = args.kind

    if not shutil.which(fc):
        print(f'Error: Fortran compiler not found: {fc}', file=sys.stderr)
        sys.exit(1)

    # Load recipe to get library name and identify common files
    config = load_recipe(args.recipe, args.project_root)
    lib_name = config.library  # e.g., "blas"

    # Determine prefix
    pmap = PREFIX_MAP[kind]
    real_pfx = pmap['D'].lower()  # q or e
    precision_lib = f'lib{real_pfx}{lib_name}.a'
    common_lib = f'lib{lib_name}_common.a'

    # Identify type-independent files via symbol classification
    symbols = scan_symbols(config.source_dir, config.language,
                           config.extensions, config.library_path)
    classification = classify_symbols(symbols)
    independent = classification.independent

    # Classify source files
    files = sorted(list(out_dir.glob('*.f')) + list(out_dir.glob('*.f90')))
    common_srcs, prec_srcs = [], []
    for f in files:
        if f.stem.upper() in independent:
            common_srcs.append(f)
        else:
            prec_srcs.append(f)

    build_dir = out_dir / '_build'
    build_dir.mkdir(exist_ok=True)

    print(f'Building {common_lib} ({len(common_srcs)} files) '
          f'+ {precision_lib} ({len(prec_srcs)} files)')

    def compile_and_archive(srcs, lib_path, prefix=''):
        objs = []
        fail = 0
        for src in srcs:
            obj = build_dir / (prefix + src.stem + '.o')
            r = subprocess.run(
                [fc] + fflags + ['-o', str(obj), str(src)],
                capture_output=True, text=True
            )
            if r.returncode == 0:
                objs.append(str(obj))
            else:
                fail += 1
                print(f'  FAIL: {src.name}')
        if fail:
            print(f'Error: {fail} file(s) failed to compile', file=sys.stderr)
            sys.exit(1)
        if objs:
            subprocess.run([ar, 'rcs', str(lib_path)] + objs, check=True)
            ranlib = shutil.which('ranlib')
            if ranlib:
                subprocess.run([ranlib, str(lib_path)], check=True)

    lib_dir = out_dir / '_lib'
    lib_dir.mkdir(exist_ok=True)

    compile_and_archive(common_srcs, lib_dir / common_lib)
    compile_and_archive(prec_srcs, lib_dir / precision_lib, prefix='p_')

    print()
    for lib in [common_lib, precision_lib]:
        p = lib_dir / lib
        if p.exists():
            # Count symbols
            r = subprocess.run(['nm', str(p)], capture_output=True, text=True)
            syms = [l for l in r.stdout.splitlines() if ' T ' in l or ' t ' in l]
            unique = set()
            for s in syms:
                name = s.split()[-1].strip('_').upper()
                if not name.startswith('LTMP'):
                    unique.add(name)
            size = p.stat().st_size
            print(f'  {lib}: {len(unique)} symbols, {size // 1024}K')

    print(f'\nLibraries in {lib_dir}/')


def cmd_run(args):
    """Run the full pipeline: migrate → verify → compile → build."""
    work_dir = args.work_dir
    output_dir = work_dir / 'output'
    output_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 60)
    print('  Step 1: Migrate')
    print('=' * 60)
    args.output_dir = output_dir
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
    print('  Step 3: Compile')
    print('=' * 60)
    args.output_dir = output_dir
    cmd_compile(args)

    print()
    print('=' * 60)
    print('  Step 4: Build')
    print('=' * 60)
    args.output_dir = output_dir
    cmd_build(args)


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
    p.set_defaults(func=cmd_migrate)

    # --- verify ---
    p = sub.add_parser('verify', help='Verify migrated output')
    p.add_argument('output_dir', type=Path)
    p.set_defaults(func=cmd_verify)

    # --- compile ---
    p = sub.add_parser('compile', help='Compile migrated files')
    p.add_argument('output_dir', type=Path)
    p.add_argument('--fc', default='gfortran', help='Fortran compiler')
    p.set_defaults(func=cmd_compile)

    # --- build ---
    p = sub.add_parser('build', help='Build static libraries')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('output_dir', type=Path)
    p.add_argument('--kind', type=int, default=16, choices=[10, 16])
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--ar', default='ar')
    p.add_argument('--project-root', type=Path, default=None)
    p.set_defaults(func=cmd_build)

    # --- run (full pipeline) ---
    p = sub.add_parser('run', help='Run full pipeline')
    p.add_argument('recipe', type=Path, help='Recipe YAML file')
    p.add_argument('work_dir', type=Path, help='Working directory')
    p.add_argument('--kind', type=int, default=16, choices=[10, 16])
    p.add_argument('--fc', default='gfortran')
    p.add_argument('--ar', default='ar')
    p.add_argument('--project-root', type=Path, default=None)
    p.set_defaults(func=cmd_run)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
