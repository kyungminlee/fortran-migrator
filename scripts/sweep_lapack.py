#!/usr/bin/env python3
"""Parallel LAPACK sweep harness for the multifloats migrator.

Migrates every eligible LAPACK source file with --target multifloats
in parallel via a multiprocessing pool, writes the output to
``out_dir``, then exits. Use ``compile_lapack.sh`` to compile the
output in parallel via ``xargs -P``.
"""

from __future__ import annotations

import shutil
import sys
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / 'src'))

from pyengine.target_mode import load_target  # noqa: E402
from pyengine.fortran_migrator import migrate_file_to_string  # noqa: E402
from pyengine.symbol_scanner import scan_symbols  # noqa: E402
from pyengine.prefix_classifier import classify_symbols  # noqa: E402
from pyengine.config import load_recipe  # noqa: E402


def _migrate_one(args):
    src_path, rename_map, target_mode = args
    try:
        return src_path, migrate_file_to_string(src_path, rename_map, target_mode, parser=None)
    except Exception as e:
        return src_path, ('ERROR', f'{type(e).__name__}: {e}')


def main() -> int:
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('/tmp/mf_lapack')
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    target_mode = load_target('multifloats')
    config = load_recipe(REPO / 'recipes' / 'lapack.yaml')
    blas_cfg = load_recipe(REPO / 'recipes' / 'blas.yaml')

    print('Scanning symbols...', flush=True)
    all_syms = list(scan_symbols(config.source_dir, 'fortran',
                                 tuple(config.extensions), None))
    all_syms += list(scan_symbols(blas_cfg.source_dir, 'fortran',
                                  tuple(blas_cfg.extensions), None))
    for d in config.extra_symbol_dirs:
        all_syms += list(scan_symbols(d, 'fortran',
                                      tuple(config.extensions), None))
    classification = classify_symbols(all_syms)
    rename_map = classification.build_rename_map(target_mode)
    # NOTE: target_mode.known_constants are handled per-file by
    # strip_known_constants_from_decls + replace_known_constants. We
    # intentionally do NOT add them to the global rename_map because
    # that would also rename the LHS aliases of LAPACK ``USE
    # LA_CONSTANTS, ONLY: zero=>dzero`` clauses, breaking the alias.
    print(f'  {len(all_syms)} symbols, {len(classification.families)} families,'
          f' {len(rename_map)} renames', flush=True)

    skip = set(config.skip_files)
    copy = set(config.copy_files)
    indep = classification.independent

    sources: list[Path] = []
    for p in sorted(config.source_dir.iterdir()):
        if p.suffix not in ('.f', '.f90', '.F90'):
            continue
        if p.stem.upper() in skip or p.stem.upper() in copy or p.stem.upper() in indep:
            continue
        sources.append(p)
    print(f'  {len(sources)} files to migrate', flush=True)

    workers = max(1, (os.cpu_count() or 4))
    print(f'  {workers} workers', flush=True)
    n_ok = n_err = 0
    errors: list[tuple[str, str]] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [
            ex.submit(_migrate_one, (p, rename_map, target_mode))
            for p in sources
        ]
        for fut in as_completed(futures):
            src_path, result = fut.result()
            if isinstance(result, tuple) and result and result[0] == 'ERROR':
                n_err += 1
                errors.append((src_path.name, result[1]))
                continue
            if result is None:
                continue
            out_name, text = result
            (out_dir / out_name).write_text(text)
            n_ok += 1

    print(f'migrated {n_ok}, errors {n_err}', flush=True)
    for name, err in errors[:10]:
        print(f'  {name}: {err}', flush=True)
    return 0 if n_err == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
