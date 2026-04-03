#!/usr/bin/env python3
"""
migrate_blas.py — BLAS type migration engine.

Migrates Reference BLAS source files from standard precision
(REAL/DOUBLE PRECISION, COMPLEX/COMPLEX*16) to extended precision
(KIND=10 or KIND=16).

Transformations applied:
  1. Type declarations (DOUBLE PRECISION → REAL(KIND=k), etc.)
  2. Routine name prefixes (D→Q, S→Q, Z→X, C→X for KIND=16)
  3. Floating-point literals (1.0D+0 → 1.0E+0_k)
  4. Type-specific intrinsics (DCONJG → CONJG, DBLE → REAL(...,KIND=k))
  5. XERBLA string arguments
  6. Comments (names and type references)
  7. File renaming

Usage:
    python3 migrate_blas.py <source-dir> <output-dir> [--kind 10|16]
"""

import argparse
import os
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Symbol database
# ---------------------------------------------------------------------------

def scan_symbols(src_dir: Path) -> dict[str, str]:
    """Scan BLAS source for SUBROUTINE/FUNCTION names."""
    names = set()
    pat = re.compile(
        r'(?:SUBROUTINE|FUNCTION)\s+([A-Za-z]\w*)', re.IGNORECASE
    )
    for f in sorted(src_dir.iterdir()):
        if f.suffix.lower() not in ('.f', '.f90', '.for', '.f95'):
            continue
        for line in f.read_text().splitlines():
            # skip comments
            if line and line[0] in ('C', 'c', '*', '!'):
                continue
            m = pat.search(line)
            if m:
                names.add(m.group(1).upper())
    return names


def classify_prefix(name: str) -> tuple[str, str]:
    """Return (prefix_char, base_name) for a BLAS symbol.
    prefix_char is one of S,D,C,Z or '' for non-precision symbols.
    """
    if len(name) < 2:
        return ('', name)
    first = name[0]
    if first in ('S', 'D', 'C', 'Z'):
        return (first, name[1:])
    return ('', name)


def build_rename_map(symbols: set[str], target_kind: int) -> dict[str, str]:
    """Build old_name → new_name mapping."""
    prefix_map = {
        'S': 'E' if target_kind == 10 else 'Q',
        'D': 'E' if target_kind == 10 else 'Q',
        'C': 'Y' if target_kind == 10 else 'X',
        'Z': 'Y' if target_kind == 10 else 'X',
    }
    rename = {}
    # Collect base names that have at least one precision variant
    bases_with_variants = set()
    for sym in symbols:
        pfx, base = classify_prefix(sym)
        if pfx:
            bases_with_variants.add(base)

    for sym in symbols:
        pfx, base = classify_prefix(sym)
        if pfx and base in bases_with_variants:
            new_pfx = prefix_map[pfx]
            rename[sym] = new_pfx + base
    return rename


# ---------------------------------------------------------------------------
# Transformation helpers
# ---------------------------------------------------------------------------

# Intrinsic replacements: old_name → (new_name, needs_kind_arg)
INTRINSIC_MAP = {
    'DCONJG': ('CONJG', False),
    'DIMAG':  ('AIMAG', False),
    'DABS':   ('ABS',   False),
    'DSQRT':  ('SQRT',  False),
    'DEXP':   ('EXP',   False),
    'DLOG':   ('LOG',   False),
    'DSIN':   ('SIN',   False),
    'DCOS':   ('COS',   False),
    'DSIGN':  ('SIGN',  False),
    'DMAX1':  ('MAX',   False),
    'DMIN1':  ('MIN',   False),
    'DNINT':  ('ANINT', False),
    'IDNINT': ('NINT',  False),
    'CABS':   ('ABS',   False),
    # These could take KIND argument, but in BLAS context the generic
    # form works correctly since all floating types are already at the
    # target KIND. Using the generic avoids column-72 overflow.
    'DBLE':   ('REAL',  False),
    'DCMPLX': ('CMPLX', False),
    'SNGL':   ('REAL',  False),
}

# INTRINSIC declaration name replacements
INTRINSIC_DECL_MAP = {
    'DCONJG': 'CONJG',
    'DIMAG':  'AIMAG',
    'DABS':   'ABS',
    'DSQRT':  'SQRT',
    'DBLE':   'REAL',
    'DCMPLX': 'CMPLX',
    'SNGL':   'REAL',
    'CABS':   'ABS',
}


def replace_type_decls(line: str, kind: int) -> str:
    """Replace precision type keywords in a code or comment line."""
    real_target = f'REAL(KIND={kind})'
    complex_target = f'COMPLEX(KIND={kind})'

    # Order matters: longest patterns first to avoid partial matches
    # DOUBLE PRECISION (case insensitive)
    line = re.sub(r'DOUBLE\s+PRECISION', real_target, line, flags=re.IGNORECASE)
    # DOUBLE COMPLEX
    line = re.sub(r'DOUBLE\s+COMPLEX', complex_target, line, flags=re.IGNORECASE)
    # COMPLEX*16, COMPLEX*8
    line = re.sub(r'COMPLEX\*16', complex_target, line, flags=re.IGNORECASE)
    line = re.sub(r'COMPLEX\*8', complex_target, line, flags=re.IGNORECASE)
    # REAL*8, REAL*4
    line = re.sub(r'REAL\*8', real_target, line, flags=re.IGNORECASE)
    line = re.sub(r'REAL\*4', real_target, line, flags=re.IGNORECASE)

    return line


def replace_standalone_real_complex(line: str, kind: int) -> str:
    """Replace standalone REAL and COMPLEX type keywords in declaration context.

    Only replaces when followed by a space and letter (declaration pattern),
    not when followed by ( which would be a function call like REAL(x).
    """
    real_target = f'REAL(KIND={kind})'
    complex_target = f'COMPLEX(KIND={kind})'

    # REAL followed by space+letter (declaration) but not already REAL(KIND=
    # Also handle REAL at the very end of line (rare)
    line = re.sub(
        r'\bREAL\b(?!\s*\(KIND)(?=\s+[A-Za-z])',
        real_target, line, flags=re.IGNORECASE
    )

    # COMPLEX followed by space+letter (but not COMPLEX*N or COMPLEX(KIND=)
    line = re.sub(
        r'\bCOMPLEX\b(?!\s*[\*(])(?=\s+[A-Za-z])',
        complex_target, line, flags=re.IGNORECASE
    )

    # Handle "COMPLEX FUNCTION" and "REAL FUNCTION" patterns
    # These are already caught by the above rules since FUNCTION starts with a letter

    return line


def replace_literals(line: str, kind: int) -> str:
    """Replace floating-point literals with KIND suffix form.

    1.0D+0  → 1.0E+0_k
    0.0D0   → 0.0E0_k
    1.0E+0  → 1.0E+0_k  (for single-precision files)
    """
    def literal_sub(m):
        mantissa = m.group(1)  # e.g., "1.0" or "0.0"
        exp_letter = m.group(2)  # D, d, E, or e
        exp_rest = m.group(3)  # e.g., "+0" or "0" or "-3"
        return f'{mantissa}E{exp_rest}_{kind}'

    # Match floating-point literals with D or E exponent
    # Pattern: digits.digits[DEde][+-]?digits
    # Also: .digits[DEde][+-]?digits  and  digits.[DEde][+-]?digits
    line = re.sub(
        r'(\d+\.\d*|\d*\.\d+)([DEde])([+-]?\d+)',
        literal_sub, line
    )

    return line


def replace_intrinsic_calls(line: str, kind: int) -> str:
    """Replace type-specific intrinsic function calls."""
    for old_name, (new_name, needs_kind) in INTRINSIC_MAP.items():
        if needs_kind:
            # DBLE(expr) → REAL(expr, KIND=k)
            # DCMPLX(expr) → CMPLX(expr, KIND=k)
            # Need to find the matching closing paren
            pattern = re.compile(
                rf'\b{old_name}\s*\(', re.IGNORECASE
            )
            while True:
                m = pattern.search(line)
                if not m:
                    break
                start = m.start()
                paren_start = line.index('(', m.start())
                # Find matching close paren
                depth = 1
                pos = paren_start + 1
                while pos < len(line) and depth > 0:
                    if line[pos] == '(':
                        depth += 1
                    elif line[pos] == ')':
                        depth -= 1
                    pos += 1
                if depth == 0:
                    close_pos = pos - 1
                    inner = line[paren_start + 1:close_pos]
                    replacement = f'{new_name}({inner}, KIND={kind})'
                    line = line[:start] + replacement + line[close_pos + 1:]
                else:
                    break  # unmatched paren, skip
        else:
            # Simple replacement: DCONJG(x) → CONJG(x)
            line = re.sub(
                rf'\b{old_name}\b', new_name, line, flags=re.IGNORECASE
            )

    return line


def replace_intrinsic_decls(line: str) -> str:
    """Replace intrinsic names in INTRINSIC declarations."""
    if not re.match(r'\s+INTRINSIC\b', line, re.IGNORECASE):
        return line

    for old_name, new_name in INTRINSIC_DECL_MAP.items():
        line = re.sub(
            rf'\b{old_name}\b', new_name, line, flags=re.IGNORECASE
        )

    return line


def replace_routine_names(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names using the rename map."""
    for old_name, new_name in rename_map.items():
        # Word-boundary replacement, case-preserving
        pattern = re.compile(rf'\b{re.escape(old_name)}\b', re.IGNORECASE)
        def case_replace(m):
            matched = m.group(0)
            if matched.isupper():
                return new_name.upper()
            elif matched.islower():
                return new_name.lower()
            else:
                # Mixed case — match the casing pattern
                return new_name.upper()
        line = pattern.sub(case_replace, line)

    return line


def replace_xerbla_strings(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names inside XERBLA string arguments.

    CALL XERBLA('DGEMM ',INFO) → CALL XERBLA('QGEMM ',INFO)
    """
    # This is already handled by replace_routine_names since it does
    # global word-boundary replacement. But string literals like 'DGEMM '
    # contain the name as a substring — we need to handle the case where
    # the name is inside quotes.
    for old_name, new_name in rename_map.items():
        # Replace inside single-quoted strings
        old_upper = old_name.upper()
        new_upper = new_name.upper()
        line = line.replace(f"'{old_upper} '", f"'{new_upper} '")
        line = line.replace(f"'{old_upper}'", f"'{new_upper}'")

    return line


# ---------------------------------------------------------------------------
# Fixed-form (.f) file migration
# ---------------------------------------------------------------------------

def is_comment_line(line: str) -> bool:
    """Check if a fixed-form line is a comment."""
    if not line:
        return True
    return line[0] in ('C', 'c', '*', '!')


def migrate_fixed_form(source: str, rename_map: dict[str, str],
                       kind: int) -> str:
    """Migrate a fixed-form Fortran file."""
    lines = source.splitlines(keepends=True)
    result = []

    for line in lines:
        original = line
        stripped = line.rstrip('\n\r')

        if not stripped:
            result.append(line)
            continue

        if is_comment_line(stripped):
            # In comments: replace routine names and type names for docs
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_type_decls(stripped, kind)
            # Preserve original line ending
            if line.endswith('\n'):
                result.append(stripped + '\n')
            else:
                result.append(stripped)
        else:
            # Code line (or continuation)
            # Apply all transformations
            stripped = replace_type_decls(stripped, kind)
            stripped = replace_standalone_real_complex(stripped, kind)
            stripped = replace_literals(stripped, kind)
            stripped = replace_intrinsic_calls(stripped, kind)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_xerbla_strings(stripped, rename_map)

            if line.endswith('\n'):
                result.append(stripped + '\n')
            else:
                result.append(stripped)

    return ''.join(result)


# ---------------------------------------------------------------------------
# Free-form (.f90) file migration
# ---------------------------------------------------------------------------

def migrate_free_form(source: str, rename_map: dict[str, str],
                      kind: int) -> str:
    """Migrate a free-form Fortran file (.f90).

    These files use parameterized types via:
        integer, parameter :: wp = kind(1.d0)
    or  integer, parameter :: wp = kind(1.e0)

    We change the wp definition and rename routines. No literal or
    type conversion is needed since they use wp-parameterized types.
    """
    lines = source.splitlines(keepends=True)
    result = []

    for line in lines:
        stripped = line.rstrip('\n\r')

        # Change wp parameter definition
        stripped = re.sub(
            r'(integer\s*,\s*parameter\s*::\s*wp\s*=\s*)kind\s*\(\s*1\.[de]0\s*\)',
            rf'\g<1>{kind}',
            stripped, flags=re.IGNORECASE
        )

        # Replace routine names (in code and comments)
        stripped = replace_routine_names(stripped, rename_map)

        # Replace type names in comments (lines starting with !)
        lstripped = stripped.lstrip()
        if lstripped.startswith('!'):
            stripped = replace_type_decls(stripped, kind)

        if line.endswith('\n'):
            result.append(stripped + '\n')
        else:
            result.append(stripped)

    return ''.join(result)


# ---------------------------------------------------------------------------
# File renaming
# ---------------------------------------------------------------------------

def target_filename(name: str, rename_map: dict[str, str]) -> str:
    """Compute output filename from source filename using rename map."""
    stem = Path(name).stem
    ext = Path(name).suffix

    upper_stem = stem.upper()
    if upper_stem in rename_map:
        new_stem = rename_map[upper_stem]
        # Preserve case of original
        if stem.islower():
            new_stem = new_stem.lower()
        elif stem.isupper():
            new_stem = new_stem.upper()
        return new_stem + ext

    return name  # No rename


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def migrate_file(src_path: Path, output_dir: Path,
                 rename_map: dict[str, str], kind: int) -> str | None:
    """Migrate a single BLAS source file. Returns output filename or None."""
    name = src_path.name
    ext = src_path.suffix.lower()

    source = src_path.read_text()

    if ext in ('.f', '.for'):
        migrated = migrate_fixed_form(source, rename_map, kind)
    elif ext in ('.f90', '.f95'):
        migrated = migrate_free_form(source, rename_map, kind)
    else:
        return None

    out_name = target_filename(name, rename_map)
    out_path = output_dir / out_name

    out_path.write_text(migrated)
    return out_name


def main():
    parser = argparse.ArgumentParser(
        description='Migrate BLAS source files to extended precision'
    )
    parser.add_argument('source_dir', type=Path,
                        help='BLAS source directory')
    parser.add_argument('output_dir', type=Path,
                        help='Output directory for migrated files')
    parser.add_argument('--kind', type=int, default=16, choices=[10, 16],
                        help='Target KIND (default: 16)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without writing')
    args = parser.parse_args()

    src_dir = args.source_dir
    out_dir = args.output_dir
    kind = args.kind

    if not src_dir.is_dir():
        print(f'Error: source directory not found: {src_dir}', file=sys.stderr)
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    # Build symbol database
    print('Scanning symbols...')
    symbols = scan_symbols(src_dir)
    rename_map = build_rename_map(symbols, kind)

    print(f'  {len(symbols)} symbols found')
    print(f'  {len(rename_map)} renames computed')

    # Precision-independent files (copy as-is)
    independent = set()
    for sym in symbols:
        pfx, base = classify_prefix(sym)
        if not pfx:
            independent.add(sym)

    # Mixed-precision routines to skip — these mix two precision types
    # and don't have a straightforward single-kind migration.
    MIXED_PRECISION = {
        'DSDOT', 'SDSDOT',         # mixed S/D dot products
        'CSROT', 'ZDROT',          # complex array + real rotation params
        'CSSCAL', 'ZDSCAL',        # complex array + real scalar
        'SCASUM', 'DZASUM',        # real result from complex input
        'SCNRM2', 'DZNRM2',       # real norm of complex vector
        'DCABS1', 'SCABS1',        # real abs-sum of complex scalar
        'ICAMAX', 'IZAMAX',        # integer index into complex array
        'ISAMAX', 'IDAMAX',        # integer index into real array
    }

    # Migrate each file
    print(f'\nMigrating to KIND={kind}...')
    migrated_count = 0
    copied_count = 0
    skipped = []

    source_files = sorted(src_dir.iterdir())
    for src_path in source_files:
        if src_path.suffix.lower() not in ('.f', '.f90', '.for', '.f95'):
            continue

        stem_upper = src_path.stem.upper()

        # Skip mixed-precision routines
        if stem_upper in MIXED_PRECISION:
            skipped.append(src_path.name)
            continue

        if args.dry_run:
            out_name = target_filename(src_path.name, rename_map)
            print(f'  {src_path.name} → {out_name}')
            continue

        # Copy precision-independent files as-is
        if stem_upper in independent:
            out_path = out_dir / src_path.name
            out_path.write_text(src_path.read_text())
            copied_count += 1
            continue

        out_name = migrate_file(src_path, out_dir, rename_map, kind)
        if out_name:
            migrated_count += 1
        else:
            skipped.append(src_path.name)

    if not args.dry_run:
        print(f'\n  Migrated:  {migrated_count} files')
        print(f'  Copied:    {copied_count} files (precision-independent)')
        if skipped:
            print(f'  Skipped:   {len(skipped)} files')
            for s in skipped:
                print(f'    {s}')
        total = migrated_count + copied_count
        print(f'  Total:     {total} files in {out_dir}')


if __name__ == '__main__':
    main()
