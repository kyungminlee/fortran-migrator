"""Shared subprocess scaffolding for the parse-tree scanners.

:mod:`flang_parser` and :mod:`gfortran_parser` differ only in which
binary they look for, which argv they build, and how they decide the
dump is usable.  Everything around that — the cached PATH search and the
missing-binary / timeout guard — is identical in both, so it lives here.
"""

from __future__ import annotations

import functools
import shutil
import subprocess
from collections.abc import Callable

# Any single dump must finish within this many seconds; a compiler that
# hangs on a pathological source must not stall a 38k-file migration.
_DUMP_TIMEOUT = 30


@functools.cache
def which_first(names: tuple[str, ...]) -> str | None:
    """Return the first of ``names`` found on PATH, or None.

    Cached — ``shutil.which`` stats every PATH entry, and the migrator
    asks once per source file when ``--parser`` is set.
    """
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return None


def run_dump(argv: list[str],
             accept: Callable[[subprocess.CompletedProcess], bool],
             ) -> str | None:
    """Run a compiler dump command and return its stdout, or None.

    ``accept`` decides whether the finished process produced a usable
    dump: the two compilers signal failure differently — flang by a
    non-zero exit status, gfortran by an empty dump (it exits 0 even
    when it only emitted warnings).
    """
    try:
        result = subprocess.run(
            argv, capture_output=True, text=True, timeout=_DUMP_TIMEOUT,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    return result.stdout if accept(result) else None
