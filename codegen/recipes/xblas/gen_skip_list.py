#!/usr/bin/env python3
"""Emit ``recipes/xblas/skip_mixed_precision.txt``.

That manifest is the file-backed spelling of ``skip_files:`` for
``recipes/xblas.yaml`` (``skip_files_manifest:``): one stem per line,
``#`` comments ignored.

XBLAS source files follow the pattern ``BLAS_<P><routine>[_<P1>][_<P2>][_x][-f2c].c``
where ``<P>`` is the output precision (s/d/c/z) and ``<P1>``/``<P2>`` are
input precisions for mixed-precision variants.

The migrator's prefix classifier groups files by their precision-letter
positions across the S/C ↔ D/Z divide. **Mixed-input variants have no
matching S↔D / C↔Z partner** in the source tree, so the classifier
cannot pair them — they would emit un-renamed under their original
stems and clash with the standard archive at link time.

These mixed-input variants are also unreachable from LAPACK (the
``*la_*_extended.f`` callers only invoke uniform-precision entry
points), so we skip them outright. This script enumerates them.

Usage::

    python recipes/xblas/gen_skip_list.py > recipes/xblas/skip_mixed_precision.txt
"""
from __future__ import annotations

import re
from pathlib import Path

# Mixed-input variant: at least one ``_<single-precision-letter>``
# segment between the routine name and the optional ``_x`` / ``-f2c``
# suffix. The leading ``BLAS_<P><routine>`` is uniform-precision; what
# follows distinguishes mixed.
#
# Pattern breakdown (lowercase, applied to stems):
#   BLAS_           literal
#   [sdcz]          output precision letter
#   [a-z][a-z0-9_]*?  routine name (non-greedy; underscores allowed
#                  for ``ge_sum_mv``)
#   (?:_[sdcz])+    one or more single-precision input segments
#   (?:_x)?         optional extra-precise suffix
#   $               end of stem
MIXED_RE = re.compile(
    r'^BLAS_[sdcz][a-z][a-z0-9_]*?(?:_[sdcz])+(?:_x)?$'
)


def collect_mixed_stems(src_root: Path) -> list[str]:
    """Return sorted list of unique mixed-input stems under src_root.

    Each .c file contributes one stem; the ``-f2c`` suffix is stripped
    so the f2c-bridge counterpart shares the stem (the migrator's
    ``skip_files`` is keyed on stems, case-insensitive).
    """
    stems: set[str] = set()
    for c_file in src_root.glob('*/BLAS_*.c'):
        stem = c_file.stem  # filename without .c
        if stem.endswith('-f2c'):
            stem = stem[:-len('-f2c')]
        if MIXED_RE.match(stem):
            stems.add(stem)
    return sorted(stems)


def main() -> None:
    src_root = Path(__file__).resolve().parents[3] / 'extern' / 'xblas-1.0.248' / 'src'
    if not src_root.is_dir():
        raise SystemExit(f'XBLAS source not found at {src_root}')

    stems = collect_mixed_stems(src_root)
    print('# Mixed-input-precision XBLAS stems skipped by the migration')
    print('# (the file-backed spelling of ``skip_files:`` — see')
    print('# recipes/xblas.yaml for why these variants are dropped).')
    print('#')
    print(f'# {len(stems)} stems, one per line; ``#`` comments and blank '
          'lines ignored,')
    print('# matched case-insensitively, without extension.')
    print('#')
    print('# Regenerate with:')
    print('#   python recipes/xblas/gen_skip_list.py '
          '> recipes/xblas/skip_mixed_precision.txt')
    for stem in stems:
        print(stem)


if __name__ == '__main__':
    main()
