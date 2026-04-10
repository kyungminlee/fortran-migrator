"""General-purpose Fortran type migration engine.

Handles fixed-form (.f) and free-form (.f90) Fortran source files.
Works for any library following the S/D/C/Z precision prefix convention
(BLAS, LAPACK, ScaLAPACK, MUMPS, etc.).

Uses a hybrid strategy (per DEVELOPER.md):
  1. Parse with Flang (via subprocess) to get a structured parse tree
  2. Extract facts: type declarations, routine names, call sites,
     literals, intrinsics, XERBLA strings
  3. Apply transformations as source-level replacements, preserving
     all original formatting, comments, and preprocessor directives

Falls back to regex-only scanning when Flang is not available.

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
from .target_mode import TargetMode


# ---------------------------------------------------------------------------
# Type declaration replacement
# ---------------------------------------------------------------------------

def replace_type_decls(line: str, target_mode: TargetMode) -> str:
    """Replace precision type keywords with target form."""
    real_target = target_mode.real_type
    complex_target = target_mode.complex_type

    # Longest patterns first to avoid partial matches
    line = re.sub(r'DOUBLE\s+PRECISION', real_target, line, flags=re.IGNORECASE)
    line = re.sub(r'DOUBLE\s+COMPLEX', complex_target, line, flags=re.IGNORECASE)
    line = re.sub(r'COMPLEX\*16', complex_target, line, flags=re.IGNORECASE)
    line = re.sub(r'COMPLEX\*8', complex_target, line, flags=re.IGNORECASE)
    line = re.sub(r'REAL\*8', real_target, line, flags=re.IGNORECASE)
    line = re.sub(r'REAL\*4', real_target, line, flags=re.IGNORECASE)

    if not target_mode.is_kind_based:
        # Check if this is a simple declaration of a known constant
        # e.g. TYPE(float64x2) ONE, ZERO
        # If so, return an empty string or comment it out
        # Matches: TYPE(name) :: VAR1, VAR2
        m = re.match(rf'^\s*{re.escape(real_target)}\s*(?:::)?\s*([A-Z0-9_,\s]+)$', line, re.IGNORECASE)
        if not m:
            m = re.match(rf'^\s*{re.escape(complex_target)}\s*(?:::)?\s*([A-Z0-9_,\s]+)$', line, re.IGNORECASE)
        
        if m:
            vars_part = m.group(1)
            vars_list = [v.strip().upper() for v in vars_part.split(',')]
            # If all variables in this line are known constants, drop the line
            if all(v in target_mode.known_constants for v in vars_list):
                return "! " + line.strip()

    return line


def replace_standalone_real_complex(line: str, target_mode: TargetMode) -> str:
    """Replace standalone REAL/COMPLEX keywords in declaration context.

    Only replaces when followed by space+letter (declaration pattern),
    not when followed by ( which would be a function call like REAL(x).
    """
    real_target = target_mode.real_type
    complex_target = target_mode.complex_type

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

def replace_literals(line: str, target_mode: TargetMode) -> str:
    """Replace floating-point literals with target form.

    KIND mode: 1.0D+0 → 1.0E+0_k, 0.0E+0 → 0.0E+0_k, 0.0 → 0.0E0_k.
    Constructor mode: 1.0D+0 → float64x2('1.0D+0') or float64x2(1.0D0).

    Bare unsuffixed literals (no D/E exponent, no ``_kind``) are also
    promoted so that complex constants like ``(0.0,0.0)`` in C sources
    converge with ``(0.0d0,0.0d0)`` in Z sources after migration.
    """
    def literal_sub(m):
        mantissa = m.group(1)
        exp_rest = m.group(3)
        if exp_rest.startswith('+'):
            exp_rest = exp_rest[1:]
            
        if target_mode.literal_mode == 'kind_suffix':
            return f'{mantissa}E{exp_rest}_{target_mode.kind_suffix}'
        else:
            orig = m.group(0)
            # Use string constructor for ALL literals in multifloats mode to avoid REAL(4) issues
            return f"{target_mode.real_constructor}('{orig}')"

    parts = re.split(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")", line)
    _FORTRAN_OP = re.compile(
        r'\.(EQ|NE|LT|GT|LE|GE|AND|OR|NOT|TRUE|FALSE|EQV|NEQV)\.',
        re.IGNORECASE,
    )
    for idx in range(0, len(parts), 2):
        seg = parts[idx]
        masked = _FORTRAN_OP.sub(
            lambda m: '\x00' + m.group(1) + '\x00', seg,
        )
        
        if target_mode.literal_mode == 'constructor':
            def _wp_sub(m):
                val = m.group(1)
                if 'D' not in val.upper() and 'E' not in val.upper():
                    if '.' in val:
                        val += 'D0'
                    else:
                        val += '.0D0'
                return f"{target_mode.real_constructor}('{val}')"
            masked = re.sub(r'(\d+\.\d*|\d*\.\d+)_wp', _wp_sub, masked, flags=re.IGNORECASE)

        masked = re.sub(
            r'(\d+\.\d*|\d*\.\d+)([DEde])([+-]?\d+)',
            literal_sub, masked,
        )
        if target_mode.literal_mode == 'kind_suffix':
            masked = re.sub(
                r'(?<![.\w])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d)',
                rf'\1E0_{target_mode.kind_suffix}', masked,
            )
        else:
            def bare_sub(m):
                val = m.group(1)
                return f"{target_mode.real_constructor}('{val}D0')"
            masked = re.sub(
                r'(?<![.\w])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d)',
                bare_sub, masked,
            )
            
        if target_mode.literal_mode == 'constructor':
            # Wrap complex constants (float64x2(...), float64x2(...)) in complex128x2(...)
            masked = re.sub(
                r'\(\s*(' + target_mode.real_constructor + r'\([^)]+\))\s*,\s*(' + target_mode.real_constructor + r'\([^)]+\))\s*\)',
                rf"{target_mode.complex_constructor}(\1,\2)",
                masked,
                flags=re.IGNORECASE
            )

        parts[idx] = re.sub(
            r'\x00([A-Za-z]+)\x00', r'.\1.', masked,
        )
    return ''.join(parts)


# ---------------------------------------------------------------------------
# Intrinsic function replacement
# ---------------------------------------------------------------------------

def replace_intrinsic_calls(line: str, target_mode: TargetMode) -> str:
    """Replace type-specific intrinsic function calls."""
    for old_name, (new_name, needs_kind) in INTRINSIC_MAP.items():
        pattern = re.compile(rf'\b{old_name}\s*\(', re.IGNORECASE)
        if needs_kind:
            search_start = 0
            while True:
                m = pattern.search(line, search_start)
                if not m:
                    break
                start = m.start()
                paren_start = line.index('(', m.start())
                depth, pos = 1, paren_start + 1
                while pos < len(line) and depth > 0:
                    if line[pos] == '(': depth += 1
                    elif line[pos] == ')': depth -= 1
                    pos += 1
                if depth == 0:
                    close_pos = pos - 1
                    inner = line[paren_start + 1:close_pos]
                    inner_stripped = inner.strip().upper()
                    old_upper = old_name.upper()
                    is_type_spec_name = old_upper in ('REAL', 'CMPLX', 'COMPLEX')
                    if (re.match(r'KIND\s*=', inner_stripped)
                            or inner_stripped.isdigit()
                            or (is_type_spec_name and re.match(r'^[A-Z_]\w*$', inner_stripped))):
                        search_start = pos
                        continue
                    
                    if target_mode.intrinsic_mode == 'add_kind':
                        depth_k = 0
                        has_top_kind = False
                        ii = 0
                        while ii < len(inner):
                            ch = inner[ii]
                            if ch == '(': depth_k += 1
                            elif ch == ')': depth_k -= 1
                            elif (depth_k == 0 and inner[ii:ii + 4].upper() == 'KIND'
                                    and re.match(r'KIND\s*=', inner[ii:], re.IGNORECASE)):
                                has_top_kind = True
                                break
                            ii += 1
                        if has_top_kind:
                            search_start = pos
                            continue
                        replacement = f'{new_name}({inner}, KIND={target_mode.kind_suffix})'
                    else:
                        # wrap_constructor mode (multifloats)
                        if new_name.upper() in ('REAL', 'AIMAG', 'CONJG'):
                            mf_name = 'MF_REAL' if new_name.upper() == 'REAL' else new_name
                            replacement = f'{mf_name}({inner})'
                        else:
                            constructor = target_mode.real_constructor if new_name.upper() == 'REAL' else target_mode.complex_constructor
                            replacement = f'{constructor}({inner})'

                    line = line[:start] + replacement + line[close_pos + 1:]
                    search_start = start + len(replacement)
                else:
                    if old_name.upper() != new_name.upper():
                        matched = line[start:m.end() - 1]
                        repl = new_name.upper() if matched.isupper() else new_name.lower()
                        line = line[:start] + repl + line[start + len(matched):]
                    break
        else:
            def _call_replace(m, _new=new_name):
                matched_name = m.group(1)
                rest = m.group(2)
                return (_new.upper() if matched_name.isupper() else _new.lower()) + rest

            line = re.sub(rf'\b({old_name})(\s*\()', _call_replace, line, flags=re.IGNORECASE)
    return line


def replace_generic_conversions(line: str, target_mode: TargetMode) -> str:
    """Add KIND (or wrap in constructor) to generic REAL() and CMPLX() calls in expression context."""
    is_fixed_cont = (len(line) >= 6 and line[:5] == '     ' and line[5] not in (' ', '0', '!', 'C', 'c', '*', '\t'))
    cont_marker = line[5] if is_fixed_cont else ''
    if is_fixed_cont:
        line = line[:5] + '\0' + line[6:]

    for name in ('REAL', 'CMPLX'):
        pattern = re.compile(rf'(?<=[=+\-*/,(.\0])\s*\b({name})\s*\(', re.IGNORECASE)
        search_start = 0
        while True:
            m = pattern.search(line, search_start)
            if not m: break
            name_start = m.start(1)
            paren_start = line.index('(', name_start)
            depth, pos = 1, paren_start + 1
            while pos < len(line) and depth > 0:
                if line[pos] == '(': depth += 1
                elif line[pos] == ')': depth -= 1
                pos += 1
            if depth != 0: break
            close_pos = pos - 1
            inner = line[paren_start + 1:close_pos]
            
            if target_mode.intrinsic_mode == 'add_kind':
                if re.search(r'\bKIND\s*=', inner, re.IGNORECASE):
                    search_start = pos
                    continue
                top_commas = 0
                d = 0
                for ch in inner:
                    if ch == '(': d += 1
                    elif ch == ')': d -= 1
                    elif ch == ',' and d == 0: top_commas += 1
                max_args = 1 if name == 'REAL' else 2
                if top_commas >= max_args:
                    search_start = pos
                    continue
                replacement = f'{name}({inner}, KIND={target_mode.kind_suffix})'
            else:
                if name.upper() == 'REAL':
                    replacement = f'MF_REAL({inner})'
                else:
                    constructor = target_mode.real_constructor if name.upper() == 'REAL' else target_mode.complex_constructor
                    replacement = f'{constructor}({inner})'

            line = line[:name_start] + replacement + line[close_pos + 1:]
            search_start = name_start + len(replacement)
    if is_fixed_cont:
        line = line[:5] + cont_marker + line[6:]
    return line


def replace_intrinsic_decls(line: str) -> str:
    """Replace intrinsic names in INTRINSIC declarations."""
    if not re.match(r'\s+INTRINSIC\b', line, re.IGNORECASE):
        return line
    for old_name, new_name in INTRINSIC_DECL_MAP.items():
        line = re.sub(rf'\b{old_name}\b', new_name, line, flags=re.IGNORECASE)

    m = re.match(r'(\s+INTRINSIC\s+)(.*)', line, re.IGNORECASE)
    if m:
        prefix, name_list = m.group(1), m.group(2)
        newline = '\n' if line.endswith('\n') else ''
        stripped = name_list.rstrip().rstrip('\n')
        trail = ''
        if stripped.endswith('&'):
            trail = ' &'
            stripped = stripped[:-1]
        elif stripped.endswith(','):
            trail = ','
            stripped = stripped[:-1]
        names = [n.strip() for n in stripped.split(',') if n.strip()]
        seen: set[str] = set()
        deduped: list[str] = []
        for n in names:
            key = n.upper()
            if key not in seen:
                seen.add(key)
                deduped.append(n)
        line = prefix + ', '.join(deduped) + trail + newline
    return line


def _dedup_intrinsic_stmts(text: str) -> str:
    """Remove duplicate names from multi-line INTRINSIC statements."""
    lines = text.split('\n')
    result: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if not re.match(r'\s+INTRINSIC\b', line, re.IGNORECASE):
            result.append(line)
            i += 1
            continue

        stmt_lines = [line]
        j = i + 1
        while j < len(lines):
            next_line = lines[j]
            is_fixed_cont = (len(next_line) > 6 and next_line[:5] == '     ' and next_line[5] not in (' ', '0', '') and next_line[0] not in ('C', 'c', '*', '!'))
            if is_fixed_cont:
                stmt_lines.append(next_line)
                j += 1
            elif stmt_lines[-1].rstrip().endswith('&'):
                stmt_lines.append(next_line)
                j += 1
            else:
                break

        m = re.match(r'(\s+INTRINSIC\s+)', stmt_lines[0], re.IGNORECASE)
        if not m:
            result.extend(stmt_lines)
            i = j
            continue

        prefix = m.group(1)
        name_part = stmt_lines[0][len(prefix):]
        for sl in stmt_lines[1:]:
            if len(sl) > 5 and sl[0] == ' ' and sl[5] not in (' ', '0', ''):
                name_part += ' ' + sl[6:]
            else:
                stripped = sl.lstrip()
                if stripped.startswith('&'): stripped = stripped[1:]
                name_part += ' ' + stripped

        name_part = name_part.replace('&', ' ')
        names = [n.strip() for n in name_part.split(',') if n.strip()]
        seen: set[str] = set()
        deduped: list[str] = []
        for n in names:
            key = n.upper()
            if key not in seen:
                seen.add(key)
                deduped.append(n)

        body = ', '.join(deduped)
        full = prefix + body
        if len(full) <= 72:
            result.append(full)
        else:
            cont_prefix = '     $                   '
            first_cap, cont_cap = 72 - len(prefix) - 1, 72 - len(cont_prefix) - 1
            chunks: list[str] = []
            cur = ''
            for name in deduped:
                addition = (', ' + name) if cur else name
                cap = first_cap if not chunks else cont_cap
                if not cur: cur = name
                elif len(cur) + len(addition) <= cap: cur += addition
                else:
                    chunks.append(cur)
                    cur = name
            if cur: chunks.append(cur)
            for ci, chunk in enumerate(chunks):
                trail = ',' if ci < len(chunks) - 1 else ''
                result.append((prefix if ci == 0 else cont_prefix) + chunk + trail)
        i = j
    return '\n'.join(result)


def replace_routine_names(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names using the rename map (case-preserving)."""
    pattern, upper_map = _get_rename_pattern(rename_map)

    def case_replace(m):
        matched = m.group(0)
        new = upper_map[matched.upper()]
        return new.upper() if matched.isupper() else (new.lower() if matched.islower() else new.upper())

    return pattern.sub(case_replace, line)

_RENAME_PATTERN_CACHE: dict[frozenset, tuple[re.Pattern, dict[str, str]]] = {}

def _get_rename_pattern(rename_map: dict[str, str]) -> tuple[re.Pattern, dict[str, str]]:
    upper_map = {k.upper(): v for k, v in rename_map.items()}
    key = frozenset(upper_map.items())
    cached = _RENAME_PATTERN_CACHE.get(key)
    if cached is not None: return cached
    names = sorted(upper_map.keys(), key=len, reverse=True)
    pattern = re.compile(r'\b(' + '|'.join(re.escape(n) for n in names) + r')\b', re.IGNORECASE) if names else re.compile(r'(?!x)x')
    _RENAME_PATTERN_CACHE[key] = (pattern, upper_map)
    return pattern, upper_map


def replace_xerbla_strings(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names inside XERBLA string arguments."""
    for old_name, new_name in rename_map.items():
        old_upper, new_upper = old_name.upper(), new_name.upper()
        line = line.replace(f"'{old_upper} '", f"'{new_upper} '")
        line = line.replace(f"'{old_upper}'", f"'{new_upper}'")
    return line


def convert_parameter_stmts(source: str, target_mode: TargetMode) -> tuple[str, list[str]]:
    """Convert floating-point PARAMETER statements to executable assignments."""
    if target_mode.is_kind_based:
        return source, []

    lines = source.splitlines(keepends=True)
    result, fp_assignments = [], []
    param_re = re.compile(r'^(\s{6,}|^\s*)PARAMETER\s*\((.*)\)\s*(!.*)?$', re.IGNORECASE)
    
    for line in lines:
        m = param_re.match(line)
        if m:
            indent, params_content, comment = m.group(1), m.group(2), m.group(3) or ''
            parts, current, depth = [], [], 0
            for char in params_content:
                if char == '(': depth += 1
                elif char == ')': depth -= 1
                if char == ',' and depth == 0:
                    parts.append(''.join(current))
                    current = []
                else: current.append(char)
            if current: parts.append(''.join(current))
            
            kept_parts = []
            for part in parts:
                if '=' in part:
                    name, val = part.split('=', 1)
                    name, val = name.strip(), val.strip()
                    is_fp = ('.' in val or 'D' in val.upper() or 'E' in val.upper() or val.upper() in target_mode.known_constants)
                    if is_fp:
                        if name.upper() in target_mode.known_constants: continue
                        fp_assignments.append(f"{indent}{name} = {val}{comment}\n")
                    else: kept_parts.append(part)
                else: kept_parts.append(part)
            
            if kept_parts: result.append(f"{indent}PARAMETER ({', '.join(kept_parts)}){comment}\n")
            else: result.append(f"{indent}! Converted to assignments: {line.strip()}\n")
        else: result.append(line)
    return "".join(result), fp_assignments


def convert_data_stmts(source: str, target_mode: TargetMode) -> tuple[str, list[str]]:
    """Convert floating-point DATA statements to executable assignments."""
    if target_mode.is_kind_based:
        return source, []

    lines = source.splitlines(keepends=True)
    result, fp_assignments = [], []
    data_re = re.compile(r'^(\s{6,}|^\s*)DATA\s+([^/]+)/\s*([^/]+)\s*/\s*(!.*)?$', re.IGNORECASE)
    
    for line in lines:
        m = data_re.match(line)
        if m:
            indent, vars_part, vals_part, comment = m.group(1), m.group(2).strip(), m.group(3).strip(), m.group(4) or ''
            if ('.' in vals_part or 'D' in vals_part.upper() or 'E' in vals_part.upper()):
                vars_list = [v.strip() for v in vars_part.split(',')]
                vals_list, current, depth = [], [], 0
                for char in vals_part:
                    if char == '(': depth += 1
                    elif char == ')': depth -= 1
                    if char == ',' and depth == 0:
                        vals_list.append(''.join(current).strip())
                        current = []
                    else: current.append(char)
                if current: vals_list.append(''.join(current).strip())
                
                if len(vars_list) == len(vals_list):
                    for v, val in zip(vars_list, vals_list):
                        if v.upper() in target_mode.known_constants: continue
                        fp_assignments.append(f"{indent}{v} = {val}{comment}\n")
                    result.append(f"{indent}! Converted to assignments: {line.strip()}\n")
                else: result.append(line)
            else: result.append(line)
        else: result.append(line)
    return "".join(result), fp_assignments


def insert_use_multifloats(source: str, target_mode: TargetMode, extra_lines: list[str] = None) -> str:
    """Insert USE multifloats statement and extra assignments after procedure headers."""
    if not target_mode.module_name and not extra_lines:
        return source

    lines = source.splitlines(keepends=True)
    result = []
    # Robust header match for SUBROUTINE, FUNCTION, PROGRAM, etc.
    proc_header_re = re.compile(
        r'^(\s{6,}|^\s*)(?:RECURSIVE\s+|PURE\s+|ELEMENTAL\s+)*'
        r'(?:(?:INTEGER|REAL|COMPLEX|LOGICAL|CHARACTER|TYPE\s*\([^)]+\)|DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX)'
        r'(?:\s*\*\s*\d+)?\s+)?'
        r'(?:PROGRAM|SUBROUTINE|FUNCTION|MODULE|BLOCK\s+DATA)\b', 
        re.IGNORECASE
    )
    
    i = 0
    while i < len(lines):
        line = lines[i]
        result.append(line)
        m = proc_header_re.match(line)
        if m:
            j = i + 1
            while j < len(lines):
                next_line = lines[j]
                if is_continuation_line(next_line):
                     result.append(next_line)
                     j += 1
                else: break
            
            if target_mode.module_name:
                indent = m.group(1)
                use_line = f"{indent if indent.strip() else '      '}USE {target_mode.module_name}\n"
                already_has = any(f"USE {target_mode.module_name}".upper() in lines[k].upper() for k in range(j, min(j+20, len(lines))))
                if not already_has: result.append(use_line)

            if extra_lines:
                k = j
                while k < len(lines):
                    l = lines[k].strip().upper()
                    if not l or l.startswith('!') or l.startswith('C') or l.startswith('*'): k += 1
                    elif any(l.startswith(p) for p in ('REAL', 'DOUBLE', 'COMPLEX', 'INTEGER', 'LOGICAL', 'CHARACTER', 'TYPE', 'USE', 'IMPLICIT', 'PARAMETER', 'DATA', 'INTRINSIC', 'EXTERNAL', 'DIMENSION', 'SAVE', 'EQUIVALENCE', 'COMMON')): k += 1
                    else: break
                for al in extra_lines: result.append(al)
                result.append("\n")
            i = j
            continue
        i += 1
    return "".join(result)


def is_comment_line(line: str) -> bool:
    return bool(line) and line[0] in ('C', 'c', '*', '!')

def is_continuation_line(line: str) -> bool:
    return len(line) > 5 and line[0:5].strip() == '' and line[5] not in (' ', '0', '')

def _build_split_mask(body: str) -> list[bool]:
    mask = [True] * len(body)
    in_string, quote_char, i = False, '', 0
    while i < len(body):
        ch = body[i]
        if in_string:
            mask[i] = False
            if ch == quote_char:
                if i + 1 < len(body) and body[i + 1] == quote_char:
                    mask[i + 1] = False
                    i += 2
                    continue
                in_string = False
        elif ch in ("'", '"'):
            in_string, quote_char, mask[i] = True, ch, False
        i += 1
    return mask

def reformat_fixed_line(line: str, cont_char: str = '+') -> str:
    if len(line) <= 72 or is_comment_line(line) or (len(line) > 6 and line[6:].lstrip().startswith('!')):
        return line
    prefix, body = line[:6] if len(line) >= 6 else line.ljust(6), line[6:]
    safe = _build_split_mask(body)
    chunks = []
    while len(body) > 66:
        split_pos = 66
        for i in range(65, max(35, 65 - 30), -1):
            if body[i] in (',', ' ') and safe[i]:
                split_pos = i + 1
                break
        else:
            for i in range(65, 0, -1):
                if safe[i]:
                    split_pos = i
                    break
        chunks.append(body[:split_pos])
        body, safe = body[split_pos:], safe[split_pos:]
    chunks.append(body)
    result_lines = [prefix + chunks[0]]
    for chunk in chunks[1:]: result_lines.append('     ' + cont_char + chunk)
    return '\n'.join(result_lines)


def migrate_fixed_form(source: str, rename_map: dict[str, str], target_mode: TargetMode) -> str:
    source, param_assignments = convert_parameter_stmts(source, target_mode)
    source, data_assignments = convert_data_stmts(source, target_mode)
    source = insert_use_multifloats(source, target_mode, extra_lines=param_assignments + data_assignments)
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        if not stripped: result.append(line); continue
        nl = '\n' if line.endswith('\n') else ''
        if is_comment_line(stripped):
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_type_decls(stripped, target_mode)
        else:
            stripped = replace_type_decls(stripped, target_mode)
            stripped = replace_standalone_real_complex(stripped, target_mode)
            stripped = replace_literals(stripped, target_mode)
            stripped = replace_intrinsic_calls(stripped, target_mode)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, target_mode)
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_xerbla_strings(stripped, rename_map)
            stripped = reformat_fixed_line(stripped)
        result.append(stripped + nl)
    
    source = ''.join(result)
    if not target_mode.is_kind_based:
        # Case-insensitive constant replacement, avoiding assignments
        for base in sorted(target_mode.known_constants.keys(), key=len, reverse=True):
            mf = target_mode.known_constants[base]
            # Replace only if not followed by = or ::
            source = re.sub(rf'(?i)\b{base}\b(?!\s*(?:=|\:\:))', mf, source)
            source = re.sub(rf'(?i)\b{mf}\b(?!\s*(?:=|\:\:))', mf, source)
            
        source = re.sub(r'(?i)\bczero\b(?!\s*(?:=|\:\:))', 'MF_ZERO', source)
        source = re.sub(r'(?i)\bsafmin\b(?!\s*(?:=|\:\:))', 'MF_SAFMIN', source)
        source = re.sub(r'(?i)\bsafmax\b(?!\s*(?:=|\:\:))', 'MF_SAFMAX', source)
        source = re.sub(r'(?i)\brtmin\b(?!\s*(?:=|\:\:))', 'MF_RTMIN', source)
        source = re.sub(r'(?i)\brtmax\b(?!\s*(?:=|\:\:))', 'MF_RTMAX', source)
        
        source = re.sub(r'! !    integer, parameter :: wp = kind\(1\.d0\)', '!    integer, parameter :: wp = kind(1.d0)', source)

    return _dedup_intrinsic_stmts(source)


_KIND_PARAM_NAMES = r'(?:wp|sp|dp)'
_KIND_PARAM_RE = re.compile(rf'(integer\s*,\s*parameter\s*::\s*{_KIND_PARAM_NAMES}\s*=\s*)(?:kind\s*\(\s*1\.[de]0\s*\)|real(?:32|64|128))', re.IGNORECASE)

def _replace_kind_parameter(line: str, target_mode: TargetMode) -> str:
    if target_mode.is_kind_based: return _KIND_PARAM_RE.sub(rf'\g<1>{target_mode.kind_suffix}', line)
    return ('! ' + line) if _KIND_PARAM_RE.search(line) else line

_ISO_USE_ONLY_RE = re.compile(r'^(?P<lead>\s*)USE\s*,\s*INTRINSIC\s*::\s*ISO_FORTRAN_ENV\s*,\s*ONLY\s*:\s*(?P<names>[^\n!]*?)\s*(?P<tail>!.*)?$', re.IGNORECASE)

def _strip_iso_fortran_env_realN(line: str) -> str:
    m = _ISO_USE_ONLY_RE.match(line)
    if not m: return line
    names = [n.strip() for n in m.group('names').split(',') if n.strip()]
    kept = [n for n in names if not re.fullmatch(r'real(?:32|64|128)', n, re.IGNORECASE)]
    if not kept: return ''
    tail = (' ' + m.group('tail')) if m.group('tail') else ''
    return f'{m.group("lead")}USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: {", ".join(kept)}{tail}'


def rewrite_la_constants_use(source: str, target_mode: TargetMode) -> str:
    const_renames = _la_constants_rename_map(target_mode)
    lines, result, in_use_stmt = source.split('\n'), [], False
    suffix = '_MF' if not target_mode.is_kind_based else '_EP'
    for line in lines:
        upper = line.upper().lstrip()
        if re.search(r'\bUSE\s+LA_XISNAN\b', upper) and f'LA_XISNAN{suffix}' not in upper:
            line = re.sub(r'(?i)\bLA_XISNAN\b', lambda m: (f'la_xisnan{suffix.lower()}' if m.group().islower() else f'LA_XISNAN{suffix}'), line)
        if re.search(r'\bUSE\s+LA_CONSTANTS\b', upper) and f'LA_CONSTANTS{suffix}' not in upper:
            in_use_stmt = True
            line = re.sub(r'(?i)\bLA_CONSTANTS\b', lambda m: (f'la_constants{suffix.lower()}' if m.group().islower() else f'LA_CONSTANTS{suffix}'), line)
        if in_use_stmt:
            line = replace_routine_names(line, const_renames)
            if not target_mode.is_kind_based:
                line = re.sub(r',\s*wp\s*=>\s*dp\s*,', ',', line, flags=re.IGNORECASE)
                line = re.sub(r',\s*wp\s*=>\s*dp\s*(?=[!&]|$)', '', line, flags=re.IGNORECASE)
                line = re.sub(r'(ONLY\s*:\s*)wp\s*=>\s*dp\s*,', r'\1', line, flags=re.IGNORECASE)
                line = re.sub(r'(ONLY\s*:\s*)wp\s*=>\s*dp\s*(?=[!&]|$)', r'\1', line, flags=re.IGNORECASE)
            if not line.rstrip().endswith('&'): in_use_stmt = False
        result.append(line)
    return '\n'.join(result)

def _la_constants_rename_map(target_mode: TargetMode) -> dict[str, str]:
    if not target_mode.is_kind_based:
        renames: dict[str, str] = {}
        for base, mf_name in target_mode.la_constants_map.items():
            u_base = base.upper()
            for p in ('S', 'D', 'C', 'Z'): renames[p + u_base] = mf_name
            renames[u_base] = mf_name
        return renames
    pmap = {10: {'S': 'E', 'D': 'E', 'C': 'Y', 'Z': 'Y'}, 16: {'S': 'Q', 'D': 'Q', 'C': 'X', 'Z': 'X'}}[target_mode.kind_suffix]
    renames: dict[str, str] = {}
    for p, names in zip(('S', 'C', 'D', 'Z'), ([ 'SP', 'SZERO', 'SHALF', 'SONE', 'STWO', 'STHREE', 'SFOUR', 'SEIGHT', 'STEN', 'SPREFIX', 'SULP', 'SEPS', 'SSAFMIN', 'SSAFMAX', 'SSMLNUM', 'SBIGNUM', 'SRTMIN', 'SRTMAX', 'STSML', 'STBIG', 'SSSML', 'SSBIG' ], ['CZERO', 'CHALF', 'CONE', 'CPREFIX'], [ 'DP', 'DZERO', 'DHALF', 'DONE', 'DTWO', 'DTHREE', 'DFOUR', 'DEIGHT', 'DTEN', 'DPREFIX', 'DULP', 'DEPS', 'DSAFMIN', 'DSAFMAX', 'DSMLNUM', 'DBIGNUM', 'DRTMIN', 'DRTMAX', 'DTSML', 'DTBIG', 'DSSML', 'DSBIG' ], ['ZZERO', 'ZHALF', 'ZONE', 'ZPREFIX'])):
        for n in names: renames[n] = pmap[p] + n[1:]
    return renames


def migrate_free_form(source: str, rename_map: dict[str, str], target_mode: TargetMode) -> str:
    source = rewrite_la_constants_use(source, target_mode)
    if not target_mode.is_kind_based:
        lines_tmp = source.splitlines()
        res_tmp = []
        in_comment_block = False
        nuke_names = {'zero', 'one', 'czero', 'safmin', 'safmax', 'rtmin', 'rtmax', 'tsml', 'tbig', 'ssml', 'sbig'}
        
        for line in lines_tmp:
            stripped = line.strip().lower()
            is_decl_start = re.match(r'^\s*(?:real|complex|integer|type|parameter).*?::', line, re.IGNORECASE) or \
                            re.match(r'^\s*parameter\s*\(', line, re.IGNORECASE)
            
            contains_nuke = False
            for n in nuke_names:
                if re.search(rf'\b{n}\b', stripped):
                    if n == 'rtmax' and 'parameter' not in stripped: continue
                    contains_nuke = True
                    break

            if not in_comment_block and is_decl_start and contains_nuke:
                res_tmp.append('! ' + line)
                if line.rstrip().endswith('&'): in_comment_block = True
            elif in_comment_block:
                res_tmp.append('! ' + line)
                if not line.rstrip().endswith('&'): in_comment_block = False
            else: res_tmp.append(line)
        source = '\n'.join(res_tmp)

    source, param_assignments = convert_parameter_stmts(source, target_mode)
    source, data_assignments = convert_data_stmts(source, target_mode)
    
    source = insert_use_multifloats(source, target_mode, extra_lines=param_assignments + data_assignments)
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        nl = '\n' if line.endswith('\n') else ''
        stripped = _replace_kind_parameter(stripped, target_mode)
        if not stripped.lstrip().startswith('!'):
            stripped = _strip_iso_fortran_env_realN(stripped)
            if not stripped: continue
            stripped = replace_intrinsic_calls(stripped, target_mode)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, target_mode)
        stripped = replace_routine_names(stripped, rename_map)
        if stripped.lstrip().startswith('!'): 
            stripped = replace_type_decls(stripped, target_mode)
        else:
            if not target_mode.is_kind_based:
                stripped = re.sub(r'REAL\(' + _KIND_PARAM_NAMES + r'\)', target_mode.real_type, stripped, flags=re.IGNORECASE)
                stripped = re.sub(r'COMPLEX\(' + _KIND_PARAM_NAMES + r'\)', target_mode.complex_type, stripped, flags=re.IGNORECASE)
        result.append(stripped + nl)
    
    source = ''.join(result)
    if not target_mode.is_kind_based:
        for base in sorted(target_mode.known_constants.keys(), key=len, reverse=True):
            mf = target_mode.known_constants[base]
            source = re.sub(rf'(?i)\b{base}\b(?!\s*(?:=|\:\:))', mf, source)
            source = re.sub(rf'(?i)\b{mf}\b(?!\s*(?:=|\:\:))', mf, source)
            
        source = re.sub(r'(?i)\bczero\b(?!\s*(?:=|\:\:))', 'MF_ZERO', source)
        source = re.sub(r'(?i)\bsafmin\b(?!\s*(?:=|\:\:))', 'MF_SAFMIN', source)
        source = re.sub(r'(?i)\bsafmax\b(?!\s*(?:=|\:\:))', 'MF_SAFMAX', source)
        source = re.sub(r'(?i)\brtmin\b(?!\s*(?:=|\:\:))', 'MF_RTMIN', source)
        source = re.sub(r'(?i)\brtmax\b(?!\s*(?:=|\:\:))', 'MF_RTMAX', source)
        source = re.sub(r'(?i)!\s*!\s*integer\s*,\s*parameter\s*::\s*wp\s*=', '!    integer, parameter :: wp =', source)
    
    return _dedup_intrinsic_stmts(source)


def target_filename(name: str, rename_map: dict[str, str]) -> str:
    stem, ext = Path(name).stem, Path(name).suffix
    if stem.upper() in rename_map:
        new = rename_map[stem.upper()]
        return (new.upper() if stem.isupper() else (new.lower() if stem.islower() else new.upper())) + ext
    return name


def migrate_file_to_string(src_path: Path, rename_map: dict[str, str], target_mode: TargetMode, parser: str | None = None, parser_cmd: str | None = None) -> tuple[str, str] | None:
    ext, source = src_path.suffix.lower(), src_path.read_text(errors='replace')
    facts = None
    if parser == 'flang':
        from .flang_parser import scan_file as flang_scan, find_flang
        cmd = parser_cmd or find_flang()
        if cmd: facts = flang_scan(src_path, cmd)
    elif parser == 'gfortran':
        from .gfortran_parser import scan_file as gfortran_scan, find_gfortran
        cmd = parser_cmd or find_gfortran()
        if cmd: facts = gfortran_scan(src_path, cmd)

    if facts is not None: migrated = _migrate_with_flang(source, ext, rename_map, target_mode, facts)
    elif ext in ('.f', '.for'): migrated = migrate_fixed_form(source, rename_map, target_mode)
    elif ext in ('.f90', '.f95', '.F90'): migrated = migrate_free_form(source, rename_map, target_mode)
    else: return None

    out_name = target_filename(src_path.name, rename_map)
    if not target_mode.is_kind_based:
        import re
        generics = {'ABS', 'REAL', 'AIMAG', 'CONJG', 'MAX', 'MIN', 'SQRT', 'EXP', 'LOG', 'LOG10', 'SIN', 'COS', 'TAN', 'SIGN', 'MOD', 'MF_REAL'}
        def clean_intrinsic(m):
            indent, funcs_str, newline = m.group(1), m.group(2), m.group(3)
            kept = [f.strip() for f in funcs_str.split(',') if f.strip() and f.strip().upper() not in generics]
            return f"{indent}INTRINSIC {', '.join(kept)}{newline}" if kept else ""
        migrated = re.sub(r'(?i)^([ \t]*)INTRINSIC\s+([A-Za-z0-9_,\s]+?)(\r?\n|$)', clean_intrinsic, migrated, flags=re.MULTILINE)
    return out_name, migrated


def migrate_file(src_path: Path, output_dir: Path, rename_map: dict[str, str], target_mode: TargetMode, parser: str | None = None, parser_cmd: str | None = None) -> str | None:
    result = migrate_file_to_string(src_path, rename_map, target_mode, parser, parser_cmd)
    if result is None: return None
    out_name, migrated = result
    (output_dir / out_name).write_text(migrated)
    return out_name


def _migrate_with_flang(source: str, ext: str, rename_map: dict[str, str], target_mode: TargetMode, facts) -> str:
    file_names = {rd.name for rd in facts.routine_defs} | {cs.name for cs in facts.call_sites} | set(facts.external_names)
    file_rename_map = {k: v for k, v in rename_map.items() if k in file_names}
    has_float_types = any(td.type_spec in ('DoublePrecision', 'Real', 'Complex') for td in facts.type_decls)
    has_real_literals = bool(facts.real_literals)
    if ext in ('.f90', '.f95', '.F90'): return _migrate_free_form_flang(source, file_rename_map, target_mode, has_float_types)
    return _migrate_fixed_form_flang(source, file_rename_map, target_mode, has_float_types, has_real_literals)


def _migrate_fixed_form_flang(source: str, rename_map: dict[str, str], target_mode: TargetMode, has_float_types: bool, has_real_literals: bool) -> str:
    source, param_assignments = convert_parameter_stmts(source, target_mode)
    source, data_assignments = convert_data_stmts(source, target_mode)
    source = insert_use_multifloats(source, target_mode, extra_lines=param_assignments + data_assignments)
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        if not stripped: result.append(line); continue
        nl = '\n' if line.endswith('\n') else ''
        if is_comment_line(stripped):
            stripped = replace_routine_names(stripped, rename_map)
            if has_float_types: stripped = replace_type_decls(stripped, target_mode)
        else:
            if has_float_types:
                stripped = replace_type_decls(stripped, target_mode)
                stripped = replace_standalone_real_complex(stripped, target_mode)
            if has_real_literals: stripped = replace_literals(stripped, target_mode)
            stripped = replace_intrinsic_calls(stripped, target_mode)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, target_mode)
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_xerbla_strings(stripped, rename_map)
            stripped = reformat_fixed_line(stripped)
        result.append(stripped + nl)
    
    source = ''.join(result)
    if not target_mode.is_kind_based:
        for base in sorted(target_mode.known_constants.keys(), key=len, reverse=True):
            mf = target_mode.known_constants[base]
            source = re.sub(rf'(?i)\b{base}\b(?!\s*(?:=|\:\:))', mf, source)
            source = re.sub(rf'(?i)\b{mf}\b(?!\s*(?:=|\:\:))', mf, source)
            
        source = re.sub(r'(?i)\bczero\b(?!\s*(?:=|\:\:))', 'MF_ZERO', source)
        source = re.sub(r'(?i)\bsafmin\b(?!\s*(?:=|\:\:))', 'MF_SAFMIN', source)
        source = re.sub(r'(?i)\bsafmax\b(?!\s*(?:=|\:\:))', 'MF_SAFMAX', source)
        source = re.sub(r'(?i)\brtmin\b(?!\s*(?:=|\:\:))', 'MF_RTMIN', source)
        source = re.sub(r'(?i)\brtmax\b(?!\s*(?:=|\:\:))', 'MF_RTMAX', source)
        
        source = re.sub(r'! !    integer, parameter :: wp = kind\(1\.d0\)', '!    integer, parameter :: wp = kind(1.d0)', source)

    return _dedup_intrinsic_stmts(source)


def _migrate_free_form_flang(source: str, rename_map: dict[str, str], target_mode: TargetMode, has_float_types: bool) -> str:
    source = rewrite_la_constants_use(source, target_mode)
    source, param_assignments = convert_parameter_stmts(source, target_mode)
    source, data_assignments = convert_data_stmts(source, target_mode)
    source = insert_use_multifloats(source, target_mode, extra_lines=param_assignments + data_assignments)
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        nl = '\n' if line.endswith('\n') else ''
        stripped = _replace_kind_parameter(stripped, target_mode)
        if not stripped.lstrip().startswith('!'):
            stripped = _strip_iso_fortran_env_realN(stripped)
            if not stripped: continue
            stripped = replace_intrinsic_calls(stripped, target_mode)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, target_mode)
        stripped = replace_routine_names(stripped, rename_map)
        if stripped.lstrip().startswith('!'): 
            stripped = replace_type_decls(stripped, target_mode)
        else:
            if not target_mode.is_kind_based:
                stripped = re.sub(r'REAL\(' + _KIND_PARAM_NAMES + r'\)', target_mode.real_type, stripped, flags=re.IGNORECASE)
                stripped = re.sub(r'COMPLEX\(' + _KIND_PARAM_NAMES + r'\)', target_mode.complex_type, stripped, flags=re.IGNORECASE)
        result.append(stripped + nl)
    
    source = ''.join(result)
    if not target_mode.is_kind_based:
        for base in sorted(target_mode.known_constants.keys(), key=len, reverse=True):
            mf = target_mode.known_constants[base]
            source = re.sub(rf'(?i)\b{base}\b(?!\s*(?:=|\:\:))', mf, source)
            source = re.sub(rf'(?i)\b{mf}\b(?!\s*(?:=|\:\:))', mf, source)
            
        source = re.sub(r'(?i)\bczero\b(?!\s*(?:=|\:\:))', 'MF_ZERO', source)
        source = re.sub(r'(?i)\bsafmin\b(?!\s*(?:=|\:\:))', 'MF_SAFMIN', source)
        source = re.sub(r'(?i)\bsafmax\b(?!\s*(?:=|\:\:))', 'MF_SAFMAX', source)
        source = re.sub(r'(?i)\brtmin\b(?!\s*(?:=|\:\:))', 'MF_RTMIN', source)
        source = re.sub(r'(?i)\brtmax\b(?!\s*(?:=|\:\:))', 'MF_RTMAX', source)
        
        source = re.sub(r'! !    integer, parameter :: wp = kind\(1\.d0\)', '!    integer, parameter :: wp = kind(1.d0)', source)
    
    return _dedup_intrinsic_stmts(source)
