"""General-purpose Fortran type migration engine.

Handles fixed-form (.f) and free-form (.f90) Fortran source files.
Works for any library following the S/D/C/Z precision prefix convention
(BLAS, LAPACK, ScaLAPACK, MUMPS, etc.).

Transformations:
  1. Type declarations (DOUBLE PRECISION → REAL(KIND=k), etc.)
  2. Standalone REAL/COMPLEX in declaration context
  3. Floating-point literals (D-exponent → E_kind suffix)
  4. Type-specific intrinsic function calls
  5. INTRINSIC declaration statements
  6. Routine name prefixes (using rename map)
  7. XERBLA string arguments
"""

import re
from pathlib import Path

from .intrinsics import INTRINSIC_MAP, INTRINSIC_DECL_MAP


# ---------------------------------------------------------------------------
# Type declaration replacement
# ---------------------------------------------------------------------------

def replace_type_decls(line: str, kind: int) -> str:
    """Replace precision type keywords with target KIND form."""
    real_target = f'REAL(KIND={kind})'
    complex_target = f'COMPLEX(KIND={kind})'

    # Longest patterns first to avoid partial matches
    line = re.sub(r'DOUBLE\s+PRECISION', real_target, line, flags=re.IGNORECASE)
    line = re.sub(r'DOUBLE\s+COMPLEX', complex_target, line, flags=re.IGNORECASE)
    line = re.sub(r'COMPLEX\*16', complex_target, line, flags=re.IGNORECASE)
    line = re.sub(r'COMPLEX\*8', complex_target, line, flags=re.IGNORECASE)
    line = re.sub(r'REAL\*8', real_target, line, flags=re.IGNORECASE)
    line = re.sub(r'REAL\*4', real_target, line, flags=re.IGNORECASE)
    return line


def replace_standalone_real_complex(line: str, kind: int) -> str:
    """Replace standalone REAL/COMPLEX keywords in declaration context.

    Only replaces when followed by space+letter (declaration pattern),
    not when followed by ( which would be a function call like REAL(x).
    """
    real_target = f'REAL(KIND={kind})'
    complex_target = f'COMPLEX(KIND={kind})'

    line = re.sub(
        r'\bREAL\b(?!\s*\(KIND)(?=\s+[A-Za-z])',
        real_target, line, flags=re.IGNORECASE
    )
    line = re.sub(
        r'\bCOMPLEX\b(?!\s*[\*(])(?=\s+[A-Za-z])',
        complex_target, line, flags=re.IGNORECASE
    )
    return line


# ---------------------------------------------------------------------------
# Literal constant replacement
# ---------------------------------------------------------------------------

def replace_literals(line: str, kind: int) -> str:
    """Replace floating-point literals with KIND suffix form.

    1.0D+0 → 1.0E+0_k, 0.0E+0 → 0.0E+0_k
    """
    def literal_sub(m):
        mantissa = m.group(1)
        exp_rest = m.group(3)
        return f'{mantissa}E{exp_rest}_{kind}'

    line = re.sub(
        r'(\d+\.\d*|\d*\.\d+)([DEde])([+-]?\d+)',
        literal_sub, line
    )
    return line


# ---------------------------------------------------------------------------
# Intrinsic function replacement
# ---------------------------------------------------------------------------

def replace_intrinsic_calls(line: str, kind: int) -> str:
    """Replace type-specific intrinsic function calls."""
    for old_name, (new_name, needs_kind) in INTRINSIC_MAP.items():
        if needs_kind:
            # Find call with matching parens and add KIND argument
            pattern = re.compile(rf'\b{old_name}\s*\(', re.IGNORECASE)
            while True:
                m = pattern.search(line)
                if not m:
                    break
                start = m.start()
                paren_start = line.index('(', m.start())
                depth, pos = 1, paren_start + 1
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
                    break
        else:
            line = re.sub(
                rf'\b{old_name}\b', new_name, line, flags=re.IGNORECASE
            )
    return line


def replace_intrinsic_decls(line: str) -> str:
    """Replace intrinsic names in INTRINSIC declarations."""
    if not re.match(r'\s+INTRINSIC\b', line, re.IGNORECASE):
        return line
    for old_name, new_name in INTRINSIC_DECL_MAP.items():
        line = re.sub(rf'\b{old_name}\b', new_name, line, flags=re.IGNORECASE)
    return line


# ---------------------------------------------------------------------------
# Routine name replacement
# ---------------------------------------------------------------------------

def replace_routine_names(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names using the rename map (case-preserving)."""
    for old_name, new_name in rename_map.items():
        pattern = re.compile(rf'\b{re.escape(old_name)}\b', re.IGNORECASE)

        def case_replace(m, _new=new_name):
            matched = m.group(0)
            if matched.isupper():
                return _new.upper()
            elif matched.islower():
                return _new.lower()
            return _new.upper()

        line = pattern.sub(case_replace, line)
    return line


def replace_xerbla_strings(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names inside XERBLA string arguments."""
    for old_name, new_name in rename_map.items():
        old_upper, new_upper = old_name.upper(), new_name.upper()
        line = line.replace(f"'{old_upper} '", f"'{new_upper} '")
        line = line.replace(f"'{old_upper}'", f"'{new_upper}'")
    return line


# ---------------------------------------------------------------------------
# Fixed-form (.f) migration
# ---------------------------------------------------------------------------

def is_comment_line(line: str) -> bool:
    """Check if a fixed-form line is a comment."""
    return bool(line) and line[0] in ('C', 'c', '*', '!')


def migrate_fixed_form(source: str, rename_map: dict[str, str],
                       kind: int) -> str:
    """Migrate a fixed-form Fortran file."""
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        if not stripped:
            result.append(line)
            continue
        nl = '\n' if line.endswith('\n') else ''
        if is_comment_line(stripped):
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_type_decls(stripped, kind)
        else:
            stripped = replace_type_decls(stripped, kind)
            stripped = replace_standalone_real_complex(stripped, kind)
            stripped = replace_literals(stripped, kind)
            stripped = replace_intrinsic_calls(stripped, kind)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_xerbla_strings(stripped, rename_map)
        result.append(stripped + nl)
    return ''.join(result)


# ---------------------------------------------------------------------------
# Free-form (.f90) migration
# ---------------------------------------------------------------------------

def migrate_free_form(source: str, rename_map: dict[str, str],
                      kind: int) -> str:
    """Migrate a free-form Fortran file (.f90).

    These files typically use parameterized types via:
        integer, parameter :: wp = kind(1.d0)
    We change the wp definition and rename routines.
    """
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        nl = '\n' if line.endswith('\n') else ''

        # Change wp parameter definition
        stripped = re.sub(
            r'(integer\s*,\s*parameter\s*::\s*wp\s*=\s*)kind\s*\(\s*1\.[de]0\s*\)',
            rf'\g<1>{kind}',
            stripped, flags=re.IGNORECASE
        )
        stripped = replace_routine_names(stripped, rename_map)

        # Replace type names in comment lines
        if stripped.lstrip().startswith('!'):
            stripped = replace_type_decls(stripped, kind)

        result.append(stripped + nl)
    return ''.join(result)


# ---------------------------------------------------------------------------
# File-level operations
# ---------------------------------------------------------------------------

def target_filename(name: str, rename_map: dict[str, str]) -> str:
    """Compute output filename from source filename using rename map."""
    stem = Path(name).stem
    ext = Path(name).suffix
    upper_stem = stem.upper()
    if upper_stem in rename_map:
        new_stem = rename_map[upper_stem]
        if stem.islower():
            new_stem = new_stem.lower()
        elif stem.isupper():
            new_stem = new_stem.upper()
        return new_stem + ext
    return name


def migrate_file(src_path: Path, output_dir: Path,
                 rename_map: dict[str, str], kind: int) -> str | None:
    """Migrate a single Fortran source file. Returns output filename."""
    ext = src_path.suffix.lower()
    source = src_path.read_text(errors='replace')

    if ext in ('.f', '.for'):
        migrated = migrate_fixed_form(source, rename_map, kind)
    elif ext in ('.f90', '.f95'):
        migrated = migrate_free_form(source, rename_map, kind)
    else:
        return None

    out_name = target_filename(src_path.name, rename_map)
    (output_dir / out_name).write_text(migrated)
    return out_name
