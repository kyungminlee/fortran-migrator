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
    """Replace precision type keywords with target form.

    In multifloats mode, also filters out variable names that are
    supplied as named constants by the multifloats module (e.g. ZERO,
    ONE, TWO). Those names become module imports and must not appear
    as locally-declared variables — otherwise gfortran complains:
    "Symbol 'mf_one' conflicts with symbol from module 'multifloats'".
    """
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
        line = _filter_known_constants_from_decl(line, target_mode)

    return line


_DECL_START_RE = re.compile(
    r'^(\s+)('
    r'DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX|'
    r'COMPLEX\s*\*\s*16|COMPLEX\s*\*\s*8|'
    r'REAL\s*\*\s*8|REAL\s*\*\s*4|'
    r'TYPE\s*\([^)]+\)|'
    # Bare REAL / COMPLEX without ``*N`` / ``(KIND=...)`` decoration.
    # The negative lookahead rules out ``REAL*N`` / ``REAL(KIND=...)``.
    # ``REAL FUNCTION FOO(...)`` is rejected later by a FUNCTION-bail
    # check on the variable list.
    r'REAL(?!\s*[*(])|COMPLEX(?!\s*[*(])'
    r')(\s*(?:,\s*[A-Za-z][\w]*\s*(?:\([^)]*\))?\s*)*::\s*|\s+)(.+)$',
    re.IGNORECASE,
)


_EQUIV_RE = re.compile(r'^\s+EQUIVALENCE\b', re.IGNORECASE | re.MULTILINE)


def _warn_on_fp_equivalence(source: str, target_mode: TargetMode) -> None:
    """Print a diagnostic if the source has an EQUIVALENCE block whose
    members would become ``TYPE(float64x2)`` after migration.

    Fortran prohibits EQUIVALENCE on non-SEQUENCE derived types, so the
    file requires manual rewriting (see Section 8 of doc/MULTIFLOATS.md).
    The migrator does not block on this — it migrates everything else
    and lets the compiler reject the EQUIVALENCE line, which gives the
    user a precise pointer to the manual fixup site.
    """
    if not _EQUIV_RE.search(source):
        return
    fp_var_pattern = re.compile(
        r'^\s+(?:DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX|'
        r'COMPLEX\s*\*\s*16|COMPLEX\s*\*\s*8|REAL\s*\*\s*[48]|'
        r'REAL(?!\s*[*(])|COMPLEX(?!\s*[*(]))\b',
        re.IGNORECASE,
    )
    fp_var_names: set[str] = set()
    for ln in source.splitlines():
        if fp_var_pattern.match(ln):
            # Extract identifiers from the variable list (rough)
            tail = re.sub(
                r'^\s+(?:DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX|'
                r'COMPLEX\s*\*\s*\d+|REAL\s*\*\s*\d+|REAL|COMPLEX)',
                '', ln, count=1, flags=re.IGNORECASE,
            )
            for m in re.finditer(r'\b([A-Za-z_]\w*)\b', tail):
                fp_var_names.add(m.group(1).upper())

    if not fp_var_names:
        return

    for ln in source.splitlines():
        m = re.match(r'^\s+EQUIVALENCE\s*(.*)$', ln, re.IGNORECASE)
        if not m:
            continue
        body = m.group(1)
        for nm in re.finditer(r'\b([A-Za-z_]\w*)\b', body):
            if nm.group(1).upper() in fp_var_names:
                import sys as _sys
                _sys.stderr.write(
                    'WARNING: EQUIVALENCE statement contains floating-point '
                    'variables that cannot be migrated to TYPE(float64x2). '
                    'Manual rewrite required (see doc/MULTIFLOATS.md §8).\n'
                    f'  -> {ln.rstrip()}\n'
                )
                return  # one warning per file is enough


def strip_known_constants_from_decls(
    source: str, target_mode: TargetMode,
) -> tuple[str, dict[str, str]]:
    """Whole-source pre-pass: drop known-constant names from multi-line type decls.

    A continuation line like ``$  DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO`` is
    invisible to per-line filtering. This pass joins continuation lines
    of a type declaration into one logical statement, removes any plain
    identifier whose uppercase form is in ``target_mode.known_constants``,
    and re-emits the declaration. Items containing parens (array specs)
    or ``=`` (initializers) are preserved verbatim. If every item is
    removed, the entire declaration is dropped.

    Returns ``(new_source, removed_renames)`` where ``removed_renames``
    maps each filtered name (uppercase) to its multifloats replacement
    (e.g. ``'ZERO' -> 'MF_ZERO'``). Callers feed this map to
    :func:`replace_known_constants` so that **only** the names that
    were actually filtered get rewritten in the body — names imported
    via ``USE LA_CONSTANTS_MF, ONLY: zero=>dzero`` are left intact
    because they were never in a local declaration.
    """
    if target_mode.is_kind_based:
        return source, {}
    known = {k.upper(): v for k, v in target_mode.known_constants.items()}
    if not known:
        return source, {}

    removed: dict[str, str] = {}
    lines = source.splitlines(keepends=True)
    out: list[str] = []
    i = 0
    while i < len(lines):
        raw = lines[i]
        # Don't touch comment lines
        if raw and raw[0] in ('C', 'c', '*', '!'):
            out.append(raw); i += 1; continue
        rstripped = raw.rstrip('\n')
        m = _DECL_START_RE.match(rstripped)
        if not m:
            out.append(raw); i += 1; continue

        indent, type_text, sep, vars_part = m.groups()

        # Bail on function-return decls: 'REAL FUNCTION FOO(...)'.
        if re.match(r'^\s*FUNCTION\b', vars_part, re.IGNORECASE):
            out.append(raw); i += 1; continue

        # Collect continuation lines (fixed-form col-6 marker, OR
        # previous logical line ends with '&').
        stmt_lines = [raw]
        j = i + 1
        prev_amp = vars_part.rstrip().endswith('&')
        while j < len(lines):
            nxt = lines[j]
            if is_continuation_line(nxt):
                stmt_lines.append(nxt); j += 1; continue
            if prev_amp:
                stmt_lines.append(nxt)
                prev_amp = nxt.rstrip('\n').rstrip().endswith('&')
                j += 1; continue
            break

        # Build the joined var-list text. Strip continuation markers
        # ('&' free-form, col-6 char fixed-form) and inline comments.
        def _strip_amp(s: str) -> str:
            s = s.rstrip()
            return s[:-1].rstrip() if s.endswith('&') else s

        full = _strip_amp(vars_part)
        for cl in stmt_lines[1:]:
            body = cl.rstrip('\n')
            if is_continuation_line(body):
                body = body[6:]
            body = body.lstrip()
            if body.startswith('&'):
                body = body[1:]
            full = full + ' ' + _strip_amp(body)

        comment = ''
        bang = full.find('!')
        if bang >= 0:
            comment = full[bang:]
            full = full[:bang]

        # Top-level comma split, respecting parentheses.
        items: list[str] = []
        cur, depth = '', 0
        for ch in full:
            if ch == '(':
                depth += 1; cur += ch
            elif ch == ')':
                depth -= 1; cur += ch
            elif ch == ',' and depth == 0:
                if cur.strip(): items.append(cur.strip())
                cur = ''
            else:
                cur += ch
        if cur.strip():
            items.append(cur.strip())

        # If any item has an '=' initializer, bail out — preserving
        # original is safer than rewriting initializer expressions.
        if any('=' in it for it in items):
            for sl in stmt_lines: out.append(sl)
            i = j; continue

        kept: list[str] = []
        for it in items:
            nm = re.match(r'^([A-Za-z_]\w*)', it)
            if nm and nm.group(1).upper() in known and it == nm.group(1):
                removed[nm.group(1).upper()] = known[nm.group(1).upper()]
                continue  # drop bare known-constant name
            kept.append(it)

        if len(kept) == len(items):
            for sl in stmt_lines: out.append(sl)
            i = j; continue

        if not kept:
            # Whole declaration removed
            i = j; continue

        body = ', '.join(kept)
        rebuilt = f'{indent}{type_text}{sep}{body}'
        if comment:
            rebuilt = rebuilt + ' ' + comment
        nl = '\n' if stmt_lines[0].endswith('\n') else ''
        out.append(rebuilt + nl)
        i = j
    return ''.join(out), removed


def _filter_known_constants_from_decl(line: str, target_mode: TargetMode) -> str:
    """Drop known-constant names from a TYPE() declaration's variable list.

    Matches `TYPE(...)` (followed by optional ``::``) followed by a
    comma-separated list of *plain* variable names — not array specs
    like ``A(LDA,*)`` or initializers, which would require a deeper
    parse. Names matching ``target_mode.known_constants`` are removed;
    if every name is removed, the entire declaration is dropped.
    """
    real_target = target_mode.real_type
    complex_target = target_mode.complex_type
    known = {k.upper() for k in target_mode.known_constants}

    type_alt = f'(?:{re.escape(real_target)}|{re.escape(complex_target)})'
    m = re.match(
        rf'^(\s*)({type_alt})(\s*(?:::)?\s*)(.+?)(\s*(?:!.*)?)$',
        line, re.IGNORECASE,
    )
    if not m:
        return line
    indent, type_text, sep, vars_part, trailer = m.groups()
    # Only handle simple comma-separated identifier lists. If any
    # entry contains parens (array spec) or '=' (initializer), bail
    # out — those need a smarter parser and the wholesale rewrite
    # would be unsafe.
    items = [v.strip() for v in vars_part.split(',')]
    if not items or any(('(' in it or '=' in it) for it in items):
        return line
    if not all(re.fullmatch(r'[A-Za-z_]\w*', it) for it in items):
        return line
    kept = [it for it in items if it.upper() not in known]
    if not kept:
        return ''  # entire declaration removed
    if len(kept) == len(items):
        return line
    return f'{indent}{type_text}{sep}{",".join(kept)}{trailer}'


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


def _all_known_constant_renames(target_mode: TargetMode) -> dict[str, str]:
    """Combined uppercase rename map for known + la_constants names."""
    out: dict[str, str] = {}
    for k, v in target_mode.known_constants.items():
        out[k.upper()] = v
    for k, v in target_mode.la_constants_map.items():
        out[k.upper()] = v
    return out


_DECL_KEYWORD_RE = re.compile(
    r'^\s*(?:TYPE\s*\(|INTEGER|REAL|COMPLEX|LOGICAL|CHARACTER|'
    r'DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX|PARAMETER\b|DATA\b|'
    r'IMPLICIT\b|DIMENSION\b|EXTERNAL\b|INTRINSIC\b|SAVE\b|'
    r'COMMON\b|EQUIVALENCE\b|USE\b)',
    re.IGNORECASE,
)


def replace_known_constants(
    line: str,
    target_mode: TargetMode,
    renames: dict[str, str] | None = None,
) -> str:
    """Substitute names in ``renames`` in the code portion of ``line`` only.

    ``renames`` is the per-file set of names that were removed from a
    local declaration by :func:`strip_known_constants_from_decls`. They
    are exactly the names that need module-import substitution. Names
    bound via a ``USE LA_CONSTANTS_MF, ONLY: zero=>dzero`` alias are
    NOT in this set and therefore are not touched here.

    When ``renames`` is ``None`` the function falls back to the union
    of ``target_mode.known_constants`` and ``target_mode.la_constants_map``;
    that legacy mode is preserved for callers that haven't yet been
    threaded through the new per-file API.

    Skips: comment lines (fixed-form column-1 ``C/c/*/!`` and bare
    ``!``-prefixed free-form), inline comment text after ``!``, and
    declaration / PARAMETER / DATA / USE statement lines (those need
    structural handling, not regex). String literals are also masked
    out so an English ``ZERO`` inside an XERBLA-style message stays
    untouched.
    """
    if target_mode.is_kind_based:
        return line
    if not line or is_comment_line(line):
        return line
    if _DECL_KEYWORD_RE.match(line):
        return line

    if renames is None:
        renames = _all_known_constant_renames(target_mode)
    if not renames:
        return line
    renames = {k.upper(): v for k, v in renames.items()}

    # Split off any inline comment so we don't substitute inside it.
    # We must respect string literals when locating the '!'.
    code_end = len(line)
    in_str, qch = False, ''
    for i, ch in enumerate(line):
        if in_str:
            if ch == qch:
                if i + 1 < len(line) and line[i + 1] == qch:
                    continue  # doubled escape
                in_str = False
        elif ch in ("'", '"'):
            in_str, qch = True, ch
        elif ch == '!':
            code_end = i
            break
    code, tail = line[:code_end], line[code_end:]

    # Mask out string literal interiors in code segment.
    parts = re.split(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")", code)
    names_alt = '|'.join(re.escape(n) for n in sorted(renames.keys(), key=len, reverse=True))
    pattern = re.compile(rf'(?<![A-Za-z0-9_])({names_alt})(?![A-Za-z0-9_])', re.IGNORECASE)

    def _sub(m):
        return renames[m.group(1).upper()]

    for idx in range(0, len(parts), 2):
        parts[idx] = pattern.sub(_sub, parts[idx])
    return ''.join(parts) + tail


def replace_xerbla_strings(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names inside XERBLA string arguments."""
    for old_name, new_name in rename_map.items():
        old_upper, new_upper = old_name.upper(), new_name.upper()
        line = line.replace(f"'{old_upper} '", f"'{new_upper} '")
        line = line.replace(f"'{old_upper}'", f"'{new_upper}'")
    return line


_FP_VALUE_RE = re.compile(
    r'^[+-]?\s*(?:'
    r'\d+\.\d*|\d*\.\d+|\d+'
    r')(?:[DdEe][+-]?\d+)?\s*$'
)


def _is_fp_value(text: str, known_constants: dict) -> bool:
    """Heuristic: does the trimmed PARAMETER/DATA value look like an
    FP literal, a known-constant reference, or an expression that
    contains one?

    Tighter than a substring scan: rejects identifiers like ``ELEMENT``
    that happen to contain ``E``. The composite-expression case splits
    on Fortran operators and checks each operand: if ANY operand is a
    known FP constant or an FP literal, the whole expression counts.
    """
    s = text.strip()
    if not s:
        return False
    if _FP_VALUE_RE.match(s):
        # Pure numeric token: FP iff it has a decimal point or exponent.
        return ('.' in s) or ('D' in s.upper()) or ('E' in s.upper())
    # Tokenize on identifiers + numeric literals; ignore operators.
    tokens = re.findall(r'[A-Za-z_]\w*|\d+\.\d*[DdEe][+-]?\d+|\d*\.\d+[DdEe]?[+-]?\d*|\d+\.\d*|\d*\.\d+', s)
    for tok in tokens:
        if tok.upper() in known_constants:
            return True
        if _FP_VALUE_RE.match(tok) and (('.' in tok) or ('D' in tok.upper()) or ('E' in tok.upper())):
            return True
    return False


def _join_continued_lines(lines: list[str], start: int) -> tuple[str, int]:
    """Join a fixed-form statement starting at ``lines[start]`` with any
    continuation lines (column-6 marker). Returns ``(joined_text, end)``
    where ``end`` is the index of the first line NOT consumed.
    """
    out = lines[start].rstrip('\n')
    j = start + 1
    while j < len(lines):
        nxt = lines[j]
        if is_continuation_line(nxt):
            out += ' ' + nxt[6:].rstrip('\n')
            j += 1
            continue
        break
    return out, j


def convert_parameter_stmts(
    source: str, target_mode: TargetMode,
) -> tuple[str, list[str], dict[str, str]]:
    """Convert floating-point PARAMETER statements to executable assignments.

    Returns ``(new_source, fp_assignments, dropped_known)`` where
    ``dropped_known`` maps each known-constant name skipped from the
    PARAMETER list to its multifloats replacement (so the caller can
    add it to the per-file rename set).

    Multi-line PARAMETER statements (fixed-form column-6 continuation)
    are joined into a single logical statement before parsing. The
    original line(s) are replaced as a unit so the line count of the
    output may differ from the input.
    """
    if target_mode.is_kind_based:
        return source, [], {}

    lines = source.splitlines(keepends=True)
    result, fp_assignments = [], []
    dropped_known: dict[str, str] = {}
    param_re = re.compile(r'^(\s{6,}|^\s*)PARAMETER\s*\((.*)\)\s*(!.*)?$', re.IGNORECASE)

    i = 0
    while i < len(lines):
        line = lines[i]
        # Try matching a single-line PARAMETER first; if not, try
        # joining continuation lines and matching the joined form.
        joined, next_i = line.rstrip('\n'), i + 1
        if param_re.match(joined) is None and re.match(r'^\s{6,}PARAMETER\b', joined, re.IGNORECASE):
            joined, next_i = _join_continued_lines(lines, i)

        m = param_re.match(joined)
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
            line_assignments: list[str] = []
            line_dropped_known: dict[str, str] = {}
            for part in parts:
                if '=' in part:
                    name, val = part.split('=', 1)
                    name, val = name.strip(), val.strip()
                    if _is_fp_value(val, target_mode.known_constants):
                        if name.upper() in target_mode.known_constants:
                            line_dropped_known[name.upper()] = target_mode.known_constants[name.upper()]
                            continue
                        line_assignments.append(f"{indent}{name} = {val}{comment}\n")
                    else: kept_parts.append(part)
                else: kept_parts.append(part)

            fp_assignments.extend(line_assignments)
            dropped_known.update(line_dropped_known)
            if kept_parts:
                result.append(f"{indent}PARAMETER ({', '.join(kept_parts)}){comment}\n")
            elif line_assignments:
                # Some FP entries became runtime assignments — leave a
                # short marker comment so reviewers can find the source.
                result.append(f"{indent}! Converted to assignments below: {joined.strip()}\n")
            # else: every entry was a known constant supplied by the
            # multifloats module — drop the line entirely (no comment).
            i = next_i
            continue
        result.append(line)
        i += 1
    return "".join(result), fp_assignments, dropped_known


def convert_data_stmts(
    source: str, target_mode: TargetMode,
) -> tuple[str, list[str], dict[str, str]]:
    """Convert floating-point DATA statements to executable assignments.

    Returns ``(new_source, fp_assignments, dropped_known)`` — see
    :func:`convert_parameter_stmts` for the meaning of the third tuple
    element.
    """
    if target_mode.is_kind_based:
        return source, [], {}

    lines = source.splitlines(keepends=True)
    result, fp_assignments = [], []
    dropped_known: dict[str, str] = {}
    data_re = re.compile(r'^(\s{6,}|^\s*)DATA\s+([^/]+)/\s*([^/]+)\s*/\s*(!.*)?$', re.IGNORECASE)

    i = 0
    while i < len(lines):
        line = lines[i]
        joined, next_i = line.rstrip('\n'), i + 1
        if data_re.match(joined) is None and re.match(r'^\s{6,}DATA\b', joined, re.IGNORECASE):
            joined, next_i = _join_continued_lines(lines, i)

        m = data_re.match(joined)
        if m:
            indent, vars_part, vals_part, comment = m.group(1), m.group(2).strip(), m.group(3).strip(), m.group(4) or ''
            # Each value is FP iff it independently matches the FP
            # heuristic; this is more discriminating than scanning the
            # whole vals_part for ``D``/``E`` substrings (which would
            # falsely match identifiers).
            tmp_vals: list[str] = []
            cur, depth = '', 0
            for ch in vals_part:
                if ch == '(': depth += 1
                elif ch == ')': depth -= 1
                if ch == ',' and depth == 0:
                    tmp_vals.append(cur.strip()); cur = ''
                else:
                    cur += ch
            if cur.strip(): tmp_vals.append(cur.strip())
            any_fp = any(_is_fp_value(v, target_mode.known_constants) for v in tmp_vals)

            if any_fp:
                vars_list = [v.strip() for v in vars_part.split(',')]
                if len(vars_list) == len(tmp_vals):
                    line_assignments: list[str] = []
                    for v, val in zip(vars_list, tmp_vals):
                        if v.upper() in target_mode.known_constants:
                            dropped_known[v.upper()] = target_mode.known_constants[v.upper()]
                            continue
                        line_assignments.append(f"{indent}{v} = {val}{comment}\n")
                    fp_assignments.extend(line_assignments)
                    if line_assignments:
                        result.append(f"{indent}! Converted to assignments below: {joined.strip()}\n")
                    # else: every name was a known constant — drop the line
                    i = next_i
                    continue
            result.append(line)
            i += 1
            continue
        result.append(line)
        i += 1
    return "".join(result), fp_assignments, dropped_known


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
                already_has = any(f"USE {target_mode.module_name}".upper() in lines[kk].upper() for kk in range(j, min(j+20, len(lines))))
                if not already_has: result.append(use_line)

            if extra_lines:
                # Walk past the declaration block (blank lines, comments,
                # and any line whose first token is a declaration keyword
                # or a fixed-form continuation of one). Insert the
                # extra assignments at the END of the declaration block,
                # i.e. just before the first executable statement —
                # otherwise the assignments would be parsed as
                # implicit-typed executables before declarations.
                k = j
                while k < len(lines):
                    raw = lines[k]
                    stripped = raw.lstrip()
                    if not stripped.strip():
                        k += 1; continue
                    if raw and raw[0] in ('C', 'c', '*', '!'):
                        k += 1; continue
                    if is_continuation_line(raw):
                        k += 1; continue
                    l = stripped.upper()
                    if any(l.startswith(p) for p in (
                        'REAL', 'DOUBLE', 'COMPLEX', 'INTEGER', 'LOGICAL',
                        'CHARACTER', 'TYPE', 'USE', 'IMPLICIT', 'PARAMETER',
                        'DATA', 'INTRINSIC', 'EXTERNAL', 'DIMENSION', 'SAVE',
                        'EQUIVALENCE', 'COMMON',
                    )):
                        k += 1; continue
                    break
                # Copy declaration block as-is, then emit assignments.
                for kk in range(j, k):
                    result.append(lines[kk])
                for al in extra_lines:
                    result.append(al)
                result.append("\n")
                i = k
            else:
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
    if not target_mode.is_kind_based:
        _warn_on_fp_equivalence(source, target_mode)
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    source, param_assignments, dropped_p = convert_parameter_stmts(source, target_mode)
    source, data_assignments, dropped_d = convert_data_stmts(source, target_mode)
    removed_known.update(dropped_p)
    removed_known.update(dropped_d)
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
            if not stripped:
                # Declaration was entirely consumed by known-constant
                # filtering — drop the line, including any newline.
                continue
            stripped = replace_standalone_real_complex(stripped, target_mode)
            stripped = replace_literals(stripped, target_mode)
            stripped = replace_intrinsic_calls(stripped, target_mode)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, target_mode)
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_xerbla_strings(stripped, rename_map)
            stripped = replace_known_constants(stripped, target_mode, renames=removed_known)
            stripped = reformat_fixed_line(stripped)
        result.append(stripped + nl)

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'! !    integer, parameter :: wp = kind\(1\.d0\)',
                        '!    integer, parameter :: wp = kind(1.d0)', source)

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
    """Rewrite ``USE LA_CONSTANTS`` clauses for the chosen target.

    KIND mode (extended precision): the LAPACK la_constants module is
    cloned to ``LA_CONSTANTS_EP`` by the migrator, so we just rename
    the module reference and rename each constant to its EP-prefixed
    equivalent (E*/Y*/Q*/X*).

    Multifloats mode: there is no ``la_constants_mf`` module — instead,
    we rewrite the import to point at the real ``multifloats`` module
    and rename each LAPACK constant (``dzero``, ``dsafmin``, ...) to its
    multifloats equivalent (``MF_ZERO``, ``MF_SAFMIN``, ...). The
    ``wp=>dp`` rename entry is dropped because ``wp`` is no longer
    meaningful once the type becomes ``TYPE(float64x2)``.
    """
    const_renames = _la_constants_rename_map(target_mode)
    lines, result, in_use_stmt = source.split('\n'), [], False
    suffix = '_MF' if not target_mode.is_kind_based else '_EP'
    target_module_upper = f'LA_CONSTANTS{suffix}'
    target_module_lower = f'la_constants{suffix.lower()}'
    target_xisnan_upper = f'LA_XISNAN{suffix}'
    target_xisnan_lower = f'la_xisnan{suffix.lower()}'

    for line in lines:
        upper = line.upper().lstrip()
        if re.search(r'\bUSE\s+LA_XISNAN\b', upper) and target_xisnan_upper not in upper:
            line = re.sub(
                r'(?i)\bLA_XISNAN\b',
                lambda m: target_xisnan_lower if m.group().islower() else target_xisnan_upper,
                line,
            )
        if re.search(r'\bUSE\s+LA_CONSTANTS\b', upper) and target_module_upper not in upper:
            in_use_stmt = True
            line = re.sub(
                r'(?i)\bLA_CONSTANTS\b',
                lambda m: target_module_lower if m.group().islower() else target_module_upper,
                line,
            )
        if in_use_stmt:
            line = replace_routine_names(line, const_renames)
            if not target_mode.is_kind_based:
                line = re.sub(r',\s*wp\s*=>\s*dp\s*,', ',', line, flags=re.IGNORECASE)
                line = re.sub(r',\s*wp\s*=>\s*dp\s*(?=[!&]|$)', '', line, flags=re.IGNORECASE)
                line = re.sub(r'(ONLY\s*:\s*)wp\s*=>\s*dp\s*,', r'\1', line, flags=re.IGNORECASE)
                line = re.sub(r'(ONLY\s*:\s*)wp\s*=>\s*dp\s*(?=[!&]|$)', r'\1', line, flags=re.IGNORECASE)
            if not line.rstrip().endswith('&'):
                in_use_stmt = False
        result.append(line)
    return '\n'.join(result)

_LA_CONSTANTS_REAL_NAMES = (
    'ZERO', 'HALF', 'ONE', 'TWO', 'THREE', 'FOUR', 'EIGHT', 'TEN',
    'PREFIX', 'ULP', 'EPS',
    'SAFMIN', 'SAFMAX', 'SMLNUM', 'BIGNUM',
    'RTMIN', 'RTMAX',
    'TSML', 'TBIG', 'SSML', 'SBIG',
)
_LA_CONSTANTS_COMPLEX_NAMES = ('ZERO', 'HALF', 'ONE', 'PREFIX')


def _la_constants_rename_map(target_mode: TargetMode) -> dict[str, str]:
    """Build a rename map for the RHS of LA_CONSTANTS USE-clause aliases.

    Maps the LAPACK la_constants names ``DZERO``, ``DSAFMIN``, ``ZZERO``,
    etc. to the equivalent names exported by the target la_constants
    auxiliary module:

      KIND=10  → ``ezero``, ``ysafmin`` (la_constants_ep)
      KIND=16  → ``qzero``, ``xsafmin`` (la_constants_ep)
      multifloats → ``wzero``, ``usafmin`` (la_constants_mf)

    Only the prefixed names are mapped — the LHS aliases ``zero``,
    ``half`` etc. are intentionally left untouched so the body of the
    routine continues to reference them through the local alias.
    """
    if target_mode.is_kind_based:
        pmap = {
            10: {'S': 'E', 'D': 'E', 'C': 'Y', 'Z': 'Y'},
            16: {'S': 'Q', 'D': 'Q', 'C': 'X', 'Z': 'X'},
        }[target_mode.kind_suffix]
    else:
        pmap = {'S': 'W', 'D': 'W', 'C': 'U', 'Z': 'U'}

    renames: dict[str, str] = {}
    for p in ('S', 'D'):
        for base in _LA_CONSTANTS_REAL_NAMES:
            renames[p + base] = pmap[p] + base
    for p in ('C', 'Z'):
        for base in _LA_CONSTANTS_COMPLEX_NAMES:
            renames[p + base] = pmap[p] + base
    return renames


def migrate_free_form(source: str, rename_map: dict[str, str], target_mode: TargetMode) -> str:
    if not target_mode.is_kind_based:
        _warn_on_fp_equivalence(source, target_mode)
    source = rewrite_la_constants_use(source, target_mode)
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    if not target_mode.is_kind_based:
        lines_tmp = source.splitlines()
        res_tmp = []
        in_comment_block = False
        # Names of free-form file-scope PARAMETER declarations that are
        # supplied as MF_* constants by the multifloats module. RTMIN and
        # RTMAX are intentionally excluded: in several LAPACK routines
        # they are local variables computed at runtime, not PARAMETERs.
        # The mapping mirrors la_constants_map but with explicit
        # multifloats target names for free-form Pattern A files (those
        # use ``USE multifloats`` directly, not ``USE LA_CONSTANTS_MF``).
        nuke_renames = {
            'zero': 'MF_ZERO', 'one': 'MF_ONE', 'czero': 'MF_ZERO',
            'safmin': 'MF_SAFMIN', 'safmax': 'MF_SAFMAX',
            'tsml': 'MF_TSML', 'tbig': 'MF_TBIG',
            'ssml': 'MF_SSML', 'sbig': 'MF_SBIG',
            # dnrm2.f90's local ``maxN = huge(0.0_wp)`` is equivalent to
            # MF_SAFMAX (which is itself defined as huge(0.0_dp) packed
            # into float64x2's high limb).
            'maxn': 'MF_SAFMAX',
        }
        nuke_names = set(nuke_renames.keys())

        for line in lines_tmp:
            stripped = line.strip().lower()
            is_decl_start = re.match(r'^\s*(?:real|complex|integer|type|parameter).*?::', line, re.IGNORECASE) or \
                            re.match(r'^\s*parameter\s*\(', line, re.IGNORECASE)

            contains_nuke = False
            matched_names: list[str] = []
            for n in nuke_names:
                if re.search(rf'\b{n}\b', stripped):
                    contains_nuke = True
                    matched_names.append(n)

            if not in_comment_block and is_decl_start and contains_nuke:
                res_tmp.append('! ' + line)
                if line.rstrip().endswith('&'): in_comment_block = True
                for n in matched_names:
                    removed_known[n.upper()] = nuke_renames[n]
            elif in_comment_block:
                res_tmp.append('! ' + line)
                if not line.rstrip().endswith('&'): in_comment_block = False
                for n in matched_names:
                    removed_known[n.upper()] = nuke_renames[n]
            else: res_tmp.append(line)
        source = '\n'.join(res_tmp)

    source, param_assignments, dropped_p = convert_parameter_stmts(source, target_mode)
    source, data_assignments, dropped_d = convert_data_stmts(source, target_mode)
    removed_known.update(dropped_p)
    removed_known.update(dropped_d)

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
            stripped = replace_known_constants(stripped, target_mode, renames=removed_known)
        result.append(stripped + nl)

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'(?i)!\s*!\s*integer\s*,\s*parameter\s*::\s*wp\s*=',
                        '!    integer, parameter :: wp =', source)

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
        # Names that the multifloats module overloads as generic
        # interfaces. INTRINSIC declarations of these names become
        # illegal once ``USE multifloats`` is in scope (gfortran:
        # "Cannot change attributes of USE-associated symbol").
        generics = {
            'ABS', 'REAL', 'AIMAG', 'CONJG', 'MAX', 'MIN', 'SQRT',
            'EXP', 'LOG', 'LOG10', 'SIN', 'COS', 'TAN', 'SIGN', 'MOD',
            'MF_REAL', 'HUGE', 'TINY', 'EPSILON', 'CEILING', 'FLOOR',
            'NINT', 'INT', 'DBLE',
            'CMPLX', 'DCMPLX', 'DCONJG', 'DIMAG',
        }
        def clean_intrinsic(m):
            indent, sep, funcs_str, newline = m.group(1), m.group(2), m.group(3), m.group(4)
            kept = [f.strip() for f in funcs_str.split(',') if f.strip() and f.strip().upper() not in generics]
            return f"{indent}INTRINSIC{sep}{', '.join(kept)}{newline}" if kept else ""
        migrated = re.sub(
            r'(?im)^([ \t]*)INTRINSIC(\s*::\s*|\s+)([A-Za-z0-9_,\s]+?)(\r?\n|$)',
            clean_intrinsic, migrated,
        )
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
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    source, param_assignments, dropped_p = convert_parameter_stmts(source, target_mode)
    source, data_assignments, dropped_d = convert_data_stmts(source, target_mode)
    removed_known.update(dropped_p)
    removed_known.update(dropped_d)
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
                if not stripped:
                    continue
                stripped = replace_standalone_real_complex(stripped, target_mode)
            if has_real_literals: stripped = replace_literals(stripped, target_mode)
            stripped = replace_intrinsic_calls(stripped, target_mode)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, target_mode)
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_xerbla_strings(stripped, rename_map)
            stripped = replace_known_constants(stripped, target_mode, renames=removed_known)
            stripped = reformat_fixed_line(stripped)
        result.append(stripped + nl)

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'! !    integer, parameter :: wp = kind\(1\.d0\)',
                        '!    integer, parameter :: wp = kind(1.d0)', source)

    return _dedup_intrinsic_stmts(source)


def _migrate_free_form_flang(source: str, rename_map: dict[str, str], target_mode: TargetMode, has_float_types: bool) -> str:
    source = rewrite_la_constants_use(source, target_mode)
    source, removed_known = strip_known_constants_from_decls(source, target_mode)
    source, param_assignments, dropped_p = convert_parameter_stmts(source, target_mode)
    source, data_assignments, dropped_d = convert_data_stmts(source, target_mode)
    removed_known.update(dropped_p)
    removed_known.update(dropped_d)
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
            stripped = replace_known_constants(stripped, target_mode, renames=removed_known)
        result.append(stripped + nl)

    source = ''.join(result)
    if not target_mode.is_kind_based:
        source = re.sub(r'(?i)!\s*!\s*integer\s*,\s*parameter\s*::\s*wp\s*=',
                        '!    integer, parameter :: wp =', source)

    return _dedup_intrinsic_stmts(source)
