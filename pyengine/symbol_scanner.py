"""Symbol scanning from source files or compiled libraries.

Extracts SUBROUTINE/FUNCTION names from Fortran source, or function
names from C source, or defined symbols from compiled archives via nm.
"""

import re
import subprocess
from pathlib import Path

# Fortran patterns
_FORTRAN_DEF_RE = re.compile(
    r'(?:SUBROUTINE|FUNCTION)\s+([A-Za-z]\w*)', re.IGNORECASE
)

# C patterns — function definitions (simplified, handles BLACS style)
_C_FUNC_RE = re.compile(
    r'^\s*(?:void|int|float|double|Int|BVOID|F_VOID_FUNC|F_INT_FUNC|F_DOUBLE_FUNC'
    r'|SCOMPLEX|DCOMPLEX|QCOMPLEX|MPI_Datatype|BLACBUFF\s*\*)'
    r'\s+(\w+)\s*\(',
    re.MULTILINE
)


def scan_fortran_source(src_dir: Path,
                        extensions: list[str] | None = None) -> set[str]:
    """Scan Fortran source files for SUBROUTINE/FUNCTION definitions."""
    if extensions is None:
        extensions = ['.f', '.f90', '.F90', '.for', '.f95']

    names: set[str] = set()
    for f in sorted(src_dir.iterdir()):
        if f.suffix.lower() not in extensions:
            continue
        for line in f.read_text(errors='replace').splitlines():
            # Skip fixed-form comments
            if line and line[0] in ('C', 'c', '*', '!'):
                continue
            # Skip free-form comments
            stripped = line.lstrip()
            if stripped.startswith('!'):
                continue
            m = _FORTRAN_DEF_RE.search(line)
            if m:
                names.add(m.group(1).upper())
    return names


def scan_c_source(src_dir: Path,
                  extensions: list[str] | None = None) -> set[str]:
    """Scan C source files for function definitions."""
    if extensions is None:
        extensions = ['.c']

    names: set[str] = set()
    for f in sorted(src_dir.iterdir()):
        if f.suffix.lower() not in extensions:
            continue
        text = f.read_text(errors='replace')
        for m in _C_FUNC_RE.finditer(text):
            names.add(m.group(1).upper())
    return names


def scan_library_archive(lib_path: Path) -> set[str]:
    """Extract defined symbols from a static or shared library using nm."""
    cmd = ['nm', '--defined-only', '-g', str(lib_path)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=30)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return set()

    nm_re = re.compile(r'\s*[0-9a-fA-F]*\s+[TtWw]\s+(\S+)')
    names: set[str] = set()
    for line in result.stdout.splitlines():
        m = nm_re.match(line)
        if m:
            sym = m.group(1)
            # Strip leading underscore (macOS) and trailing underscore (Fortran)
            if sym.startswith('_'):
                sym = sym[1:]
            if sym.endswith('_'):
                sym = sym[:-1]
            if sym:
                names.add(sym.upper())
    return names


def scan_symbols(source_dir: Path, language: str,
                 extensions: list[str] | None = None,
                 library_path: Path | None = None) -> set[str]:
    """Scan for symbols using the appropriate method.

    If library_path is provided, uses nm. Otherwise scans source.
    """
    if library_path and library_path.exists():
        return scan_library_archive(library_path)

    if language == 'fortran':
        return scan_fortran_source(source_dir, extensions)
    elif language == 'c':
        return scan_c_source(source_dir, extensions)
    else:
        raise ValueError(f'Unsupported language: {language}')
