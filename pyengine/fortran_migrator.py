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

    1.0D+0 → 1.0E+0_k, 0.0E+0 → 0.0E+0_k, 0.0 → 0.0E0_k.

    Bare unsuffixed literals (no D/E exponent, no ``_kind``) are also
    promoted so that complex constants like ``(0.0,0.0)`` in C sources
    converge with ``(0.0d0,0.0d0)`` in Z sources after migration to
    KIND=16.
    """
    def literal_sub(m):
        mantissa = m.group(1)
        exp_rest = m.group(3)
        # Drop a leading ``+`` (canonical form: ``E0``, not ``E+0``)
        # so C sources like ``0.0E+0`` converge with Z sources like
        # ``0.0d0`` after migration.
        if exp_rest.startswith('+'):
            exp_rest = exp_rest[1:]
        return f'{mantissa}E{exp_rest}_{kind}'

    # Split off string literals so FORMAT specs like ``(E12.4)`` are
    # never rewritten; rejoin after processing each code segment.
    parts = re.split(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")", line)
    # Fortran logical operators ``.EQ.``/``.AND.``/... end in ``.``;
    # their trailing dot must not fool the bare-literal regex into
    # matching ``.EQ.0.0`` as ``.EQ`` + ``.0`` + ``.0``. Neutralize
    # the operator bodies with a placeholder that has no ``.``, then
    # restore once the literal rewrite is done.
    _FORTRAN_OP = re.compile(
        r'\.(EQ|NE|LT|GT|LE|GE|AND|OR|NOT|TRUE|FALSE|EQV|NEQV)\.',
        re.IGNORECASE,
    )
    for idx in range(0, len(parts), 2):
        seg = parts[idx]
        masked = _FORTRAN_OP.sub(
            lambda m: '\x00' + m.group(1) + '\x00', seg,
        )
        masked = re.sub(
            r'(\d+\.\d*|\d*\.\d+)([DEde])([+-]?\d+)',
            literal_sub, masked,
        )
        # Bare float (no exponent, no kind suffix). Lookbehind avoids
        # identifiers/digits; ``.`` cannot precede a literal here
        # because logical operators are already masked. Lookahead
        # stops us re-processing literals already rewritten above
        # (``0.0E0_16``) or touching identifier / kind suffixes.
        masked = re.sub(
            r'(?<![.\w])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d)',
            rf'\1E0_{kind}', masked,
        )
        parts[idx] = re.sub(
            r'\x00([A-Za-z]+)\x00', r'.\1.', masked,
        )
    return ''.join(parts)


# ---------------------------------------------------------------------------
# Intrinsic function replacement
# ---------------------------------------------------------------------------

def replace_intrinsic_calls(line: str, kind: int) -> str:
    """Replace type-specific intrinsic function calls.

    Only replaces when followed by '(' (function call syntax) to avoid
    renaming variables that happen to share names with intrinsics
    (e.g., DMIN1 as a local variable vs. DMIN1() intrinsic).
    """
    for old_name, (new_name, needs_kind) in INTRINSIC_MAP.items():
        pattern = re.compile(rf'\b{old_name}\s*\(', re.IGNORECASE)
        if needs_kind:
            # Find call with matching parens and add KIND argument
            search_start = 0
            while True:
                m = pattern.search(line, search_start)
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
                    # Skip type declarations: REAL(KIND=16), REAL(4),
                    # REAL(wp), COMPLEX(dp), etc.
                    # The single-identifier check only applies when the
                    # old name can itself be a Fortran type specifier.
                    # Intrinsics like DBLE, SNGL, DCMPLX are never type
                    # specifiers, so DBLE(ALPHA) is always a function call.
                    inner_stripped = inner.strip().upper()
                    old_upper = old_name.upper()
                    is_type_spec_name = old_upper in ('REAL', 'CMPLX',
                                                       'COMPLEX')
                    if (re.match(r'KIND\s*=', inner_stripped)
                            or inner_stripped.isdigit()
                            or (is_type_spec_name
                                and re.match(r'^[A-Z_]\w*$',
                                             inner_stripped))):
                        search_start = pos
                        continue
                    # Skip if this call already carries a ``KIND=``
                    # argument *at top level* (meaning we already added
                    # it to this very call). Nested ``KIND=`` inside a
                    # sub-call — e.g. ``DCMPLX(DA*REAL(X, KIND=16),…)``
                    # after the inner ``DBLE``→``REAL`` rewrite — must
                    # not block the outer intrinsic's own rename.
                    depth_k = 0
                    has_top_kind = False
                    ii = 0
                    while ii < len(inner):
                        ch = inner[ii]
                        if ch == '(':
                            depth_k += 1
                        elif ch == ')':
                            depth_k -= 1
                        elif (depth_k == 0
                                and inner[ii:ii + 4].upper() == 'KIND'
                                and re.match(r'KIND\s*=',
                                             inner[ii:],
                                             re.IGNORECASE)):
                            has_top_kind = True
                            break
                        ii += 1
                    if has_top_kind:
                        search_start = pos
                        continue
                    replacement = f'{new_name}({inner}, KIND={kind})'
                    line = line[:start] + replacement + line[close_pos + 1:]
                    search_start = start + len(replacement)
                else:
                    # Call spans continuation lines — just replace the name
                    if old_name.upper() != new_name.upper():
                        matched = line[start:m.end() - 1]
                        if matched.isupper():
                            repl = new_name.upper()
                        elif matched.islower():
                            repl = new_name.lower()
                        else:
                            repl = new_name.upper()
                        line = line[:start] + repl + line[start + len(matched):]
                    break
        else:
            # Replace only the name part before '(', not bare occurrences
            def _call_replace(m, _new=new_name):
                matched_name = m.group(1)
                rest = m.group(2)  # whitespace + '('
                if matched_name.isupper():
                    return _new.upper() + rest
                elif matched_name.islower():
                    return _new.lower() + rest
                return _new.upper() + rest

            line = re.sub(
                rf'\b({old_name})(\s*\()', _call_replace, line,
                flags=re.IGNORECASE
            )
    return line


def replace_generic_conversions(line: str, kind: int) -> str:
    """Add KIND to generic REAL() and CMPLX() calls in expression context.

    Distinguishes function calls from type declarations by checking
    what precedes the call: operators/commas/parens indicate expression
    context, while line-start indicates a type declaration.
    """
    # Pattern: REAL or CMPLX preceded by expression-context character.
    # The '.' handles Fortran relational operators (.EQ./.NE./.NOT./etc.)
    # that precede REAL() as a type-conversion intrinsic call.
    #
    # Fixed-form continuation lines begin with a non-blank marker in
    # column 6 followed by whitespace and then the continuation of the
    # previous statement's expression. For these the whitespace after
    # column 6 is expression context, so ``REAL( TJJS )`` at the start
    # of a continuation body is a type-conversion call, not a type
    # declaration. The lookbehind below includes ``\0`` (NUL) as an
    # expression-context marker; if the line is a fixed-form
    # continuation we temporarily replace its column-6 marker with NUL
    # so the lookbehind matches, then restore it before returning.
    # Every edit we make is strictly past column 6, so this swap is
    # invariant with respect to the substitutions below.
    is_fixed_cont = (
        len(line) >= 6 and line[:5] == '     '
        and line[5] not in (' ', '0', '!', 'C', 'c', '*', '\t')
    )
    cont_marker = line[5] if is_fixed_cont else ''
    if is_fixed_cont:
        line = line[:5] + '\0' + line[6:]

    for name in ('REAL', 'CMPLX'):
        pattern = re.compile(
            rf'(?<=[=+\-*/,(.\0])\s*\b({name})\s*\(', re.IGNORECASE
        )
        search_start = 0
        while True:
            m = pattern.search(line, search_start)
            if not m:
                break
            # Find the name start and paren
            name_start = m.start(1)
            paren_start = line.index('(', name_start)
            depth, pos = 1, paren_start + 1
            while pos < len(line) and depth > 0:
                if line[pos] == '(':
                    depth += 1
                elif line[pos] == ')':
                    depth -= 1
                pos += 1
            if depth != 0:
                break  # multi-line call, skip
            close_pos = pos - 1
            inner = line[paren_start + 1:close_pos]
            if re.search(r'\bKIND\s*=', inner, re.IGNORECASE):
                search_start = pos
                continue
            # Count top-level commas (not inside nested parens)
            top_commas = 0
            d = 0
            for ch in inner:
                if ch == '(':
                    d += 1
                elif ch == ')':
                    d -= 1
                elif ch == ',' and d == 0:
                    top_commas += 1
            # REAL(x, kind) already has kind as 2nd arg; skip
            # CMPLX(x, y, kind) already has kind as 3rd arg; skip
            max_args = 1 if name == 'REAL' else 2
            if top_commas >= max_args:
                search_start = pos
                continue
            replacement = f'{name}({inner}, KIND={kind})'
            line = line[:name_start] + replacement + line[close_pos + 1:]
            search_start = name_start + len(replacement)
    if is_fixed_cont:
        line = line[:5] + cont_marker + line[6:]
    return line


def replace_intrinsic_decls(line: str) -> str:
    """Replace intrinsic names in INTRINSIC declarations.

    After substitution, removes duplicate names that can arise when a
    type-specific intrinsic (e.g. DBLE) maps to a generic (REAL) that
    already appears in the same declaration.
    """
    if not re.match(r'\s+INTRINSIC\b', line, re.IGNORECASE):
        return line
    for old_name, new_name in INTRINSIC_DECL_MAP.items():
        line = re.sub(rf'\b{old_name}\b', new_name, line, flags=re.IGNORECASE)

    # Deduplicate: parse out the name list, remove duplicates, reassemble.
    m = re.match(r'(\s+INTRINSIC\s+)(.*)', line, re.IGNORECASE)
    if m:
        prefix, name_list = m.group(1), m.group(2)
        newline = '\n' if line.endswith('\n') else ''
        stripped = name_list.rstrip().rstrip('\n')
        # Detect trailing continuation marker (& for free-form,
        # trailing comma for fixed-form continuation)
        trail = ''
        if stripped.endswith('&'):
            trail = ' &'
            stripped = stripped[:-1]
        elif stripped.endswith(','):
            # Trailing comma indicates a continuation line follows
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


# ---------------------------------------------------------------------------
# Routine name replacement
# ---------------------------------------------------------------------------

_RENAME_PATTERN_CACHE: dict[
    frozenset, tuple[re.Pattern, dict[str, str]]
] = {}


def _get_rename_pattern(
    rename_map: dict[str, str],
) -> tuple[re.Pattern, dict[str, str]]:
    """Build (and cache) a single alternation regex over all rename keys.

    A single pattern is 10–100× faster than iterating every key for
    every line — dominant cost on LAPACK (2000+ renames × many
    thousand lines).

    Cache keyed by the frozen set of ``(upper_key, value)`` pairs so
    two filtered rename maps with the same contents share a compiled
    pattern. Keying by ``id()`` is unsafe because Python reuses
    addresses after a dict is freed, which caused stale-pattern
    collisions when ``_migrate_with_flang`` built a fresh
    ``file_rename_map`` per file."""
    upper_map = {k.upper(): v for k, v in rename_map.items()}
    key = frozenset(upper_map.items())
    cached = _RENAME_PATTERN_CACHE.get(key)
    if cached is not None:
        return cached
    # Longest first so a key that's a prefix of another never steals
    # the match from its longer sibling.
    names = sorted(upper_map.keys(), key=len, reverse=True)
    if not names:
        pattern = re.compile(r'(?!x)x')  # never matches
    else:
        pattern = re.compile(
            r'\b(' + '|'.join(re.escape(n) for n in names) + r')\b',
            re.IGNORECASE,
        )
    _RENAME_PATTERN_CACHE[key] = (pattern, upper_map)
    return pattern, upper_map


def replace_routine_names(line: str, rename_map: dict[str, str]) -> str:
    """Replace routine names using the rename map (case-preserving)."""
    pattern, upper_map = _get_rename_pattern(rename_map)

    def case_replace(m):
        matched = m.group(0)
        new = upper_map[matched.upper()]
        if matched.isupper():
            return new.upper()
        elif matched.islower():
            return new.lower()
        return new.upper()

    return pattern.sub(case_replace, line)


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


def is_continuation_line(line: str) -> bool:
    """Check if a fixed-form line is a continuation."""
    return len(line) > 5 and line[0:5].strip() == '' and line[5] not in (' ', '0', '')


def reformat_fixed_line(line: str, cont_char: str = '+') -> str:
    """Split a fixed-form code line that exceeds column 72.

    Splits the statement portion into chunks that fit within columns 7-72
    (66 chars per line), using continuation lines with the given character
    in column 6.
    """
    if len(line) <= 72:
        return line
    if is_comment_line(line):
        return line  # Don't split comments
    # Don't split lines where everything past column 6 is a ! comment
    if len(line) > 6 and line[6:].lstrip().startswith('!'):
        return line

    # Extract label/continuation prefix (columns 1-6) and statement body
    prefix = line[:6] if len(line) >= 6 else line.ljust(6)
    body = line[6:]

    # Split body into 66-char chunks
    chunks = []
    while len(body) > 66:
        # Try to break at a comma or space
        split_pos = 66
        for i in range(65, max(19, 65 - 30), -1):
            if body[i] in (',', ' '):
                split_pos = i + 1
                break
        chunks.append(body[:split_pos])
        body = body[split_pos:]
    chunks.append(body)

    # Reassemble: first chunk with original prefix, rest with continuation
    result_lines = [prefix + chunks[0]]
    cont_prefix = '     ' + cont_char
    for chunk in chunks[1:]:
        result_lines.append(cont_prefix + chunk)

    return '\n'.join(result_lines)


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
            stripped = replace_generic_conversions(stripped, kind)
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_xerbla_strings(stripped, rename_map)
            # Handle column 72 overflow
            stripped = reformat_fixed_line(stripped)
        result.append(stripped + nl)
    return ''.join(result)


# ---------------------------------------------------------------------------
# Free-form (.f90) migration
# ---------------------------------------------------------------------------

# Regex for kind parameter definitions:
#   integer, parameter :: wp = kind(1.d0)
#   integer, parameter :: sp = kind(1.e0)
#   integer, parameter :: dp = real64
_KIND_PARAM_NAMES = r'(?:wp|sp|dp)'
_KIND_PARAM_VALUES = r'(?:kind\s*\(\s*1\.[de]0\s*\)|real(?:32|64|128))'
_KIND_PARAM_RE = re.compile(
    rf'(integer\s*,\s*parameter\s*::\s*{_KIND_PARAM_NAMES}\s*=\s*){_KIND_PARAM_VALUES}',
    re.IGNORECASE,
)


def _replace_kind_parameter(line: str, kind: int) -> str:
    """Replace working-precision parameter definitions."""
    return _KIND_PARAM_RE.sub(rf'\g<1>{kind}', line)


# ``USE, INTRINSIC :: iso_fortran_env, ONLY: real{32,64,128}`` survives
# migration even though the matching ``INTEGER, PARAMETER :: WP = real64``
# is rewritten to ``WP = {kind}`` — so the import becomes dead code and
# the precision-suffix carries the original (non-converged) realN name.
# Strip the realN entry from the ONLY list, and delete the whole
# statement when realN was its sole import.
_ISO_USE_ONLY_RE = re.compile(
    r'^(?P<lead>\s*)USE\s*,\s*INTRINSIC\s*::\s*ISO_FORTRAN_ENV'
    r'\s*,\s*ONLY\s*:\s*(?P<names>[^\n!]*?)\s*(?P<tail>!.*)?$',
    re.IGNORECASE,
)


def _strip_iso_fortran_env_realN(line: str) -> str:
    """Drop ``real{32,64,128}`` from an ``iso_fortran_env`` ONLY list.

    After ``_replace_kind_parameter`` rewrites ``WP = real64`` to
    ``WP = 16`` the import is unused. Trim it so co-family halves
    whose originals referenced different realN names converge.
    """
    m = _ISO_USE_ONLY_RE.match(line)
    if not m:
        return line
    names = [n.strip() for n in m.group('names').split(',') if n.strip()]
    kept = [n for n in names
            if not re.fullmatch(r'real(?:32|64|128)', n, re.IGNORECASE)]
    if not kept:
        return ''   # drop the whole USE statement
    tail = m.group('tail') or ''
    tail = (' ' + tail) if tail else ''
    return f'{m.group("lead")}USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: {", ".join(kept)}{tail}'


# ---------------------------------------------------------------------------
# USE LA_CONSTANTS → USE LA_CONSTANTS_EP rewriting
# ---------------------------------------------------------------------------

# All parameter names defined in la_constants.f90, grouped by prefix.
_LA_CONST_SP = [
    'SP', 'SZERO', 'SHALF', 'SONE', 'STWO', 'STHREE', 'SFOUR',
    'SEIGHT', 'STEN', 'SPREFIX',
    'SULP', 'SEPS', 'SSAFMIN', 'SSAFMAX', 'SSMLNUM', 'SBIGNUM',
    'SRTMIN', 'SRTMAX', 'STSML', 'STBIG', 'SSSML', 'SSBIG',
]
_LA_CONST_CP = ['CZERO', 'CHALF', 'CONE', 'CPREFIX']
_LA_CONST_DP = [
    'DP', 'DZERO', 'DHALF', 'DONE', 'DTWO', 'DTHREE', 'DFOUR',
    'DEIGHT', 'DTEN', 'DPREFIX',
    'DULP', 'DEPS', 'DSAFMIN', 'DSAFMAX', 'DSMLNUM', 'DBIGNUM',
    'DRTMIN', 'DRTMAX', 'DTSML', 'DTBIG', 'DSSML', 'DSBIG',
]
_LA_CONST_ZP = ['ZZERO', 'ZHALF', 'ZONE', 'ZPREFIX']


def _la_constants_rename_map(kind: int) -> dict[str, str]:
    """Build a rename map for LA_CONSTANTS symbols at a given KIND.

    LA_CONSTANTS exposes precision-tagged constants (SZERO/DZERO/CZERO/
    ZZERO and friends) that all migrate to the extended-precision name.
    This is a special case — unlike routine renames, all four prefix
    variants must be remapped, so we use a local S/C/D/Z → Q/X (or
    E/Y) map independent of the generic PREFIX_MAP.
    """
    _LA_FULL_MAP = {
        10: {'S': 'E', 'D': 'E', 'C': 'Y', 'Z': 'Y'},
        16: {'S': 'Q', 'D': 'Q', 'C': 'X', 'Z': 'X'},
    }
    pmap = _LA_FULL_MAP[kind]
    renames: dict[str, str] = {}
    for name in _LA_CONST_SP:
        renames[name] = pmap['S'] + name[1:]
    for name in _LA_CONST_CP:
        renames[name] = pmap['C'] + name[1:]
    for name in _LA_CONST_DP:
        renames[name] = pmap['D'] + name[1:]
    for name in _LA_CONST_ZP:
        renames[name] = pmap['Z'] + name[1:]
    return renames


def rewrite_la_constants_use(source: str, kind: int) -> str:
    """Rewrite USE LA_CONSTANTS → LA_CONSTANTS_EP and USE LA_XISNAN → LA_XISNAN_EP.

    Handles multi-line USE statements (free-form & continuation).
    For LA_CONSTANTS: also renames imported constant symbols to their
    extended/quad-precision counterparts.
    """
    const_renames = _la_constants_rename_map(kind)
    lines = source.split('\n')
    result: list[str] = []
    in_use_stmt = False  # tracks multi-line LA_CONSTANTS USE

    for line in lines:
        upper = line.upper().lstrip()

        # --- USE LA_XISNAN → USE LA_XISNAN_EP (always single-line) ---
        if re.search(r'\bUSE\s+LA_XISNAN\b', upper) and 'LA_XISNAN_EP' not in upper:
            line = re.sub(
                r'(?i)\bLA_XISNAN\b',
                lambda m: ('la_xisnan_ep' if m.group().islower()
                           else 'LA_XISNAN_EP'),
                line,
            )

        # --- USE LA_CONSTANTS → USE LA_CONSTANTS_EP + rename symbols ---
        if re.search(r'\bUSE\s+LA_CONSTANTS\b', upper) and 'LA_CONSTANTS_EP' not in upper:
            in_use_stmt = True
            line = re.sub(
                r'(?i)\bLA_CONSTANTS\b',
                lambda m: ('la_constants_ep' if m.group().islower()
                           else 'LA_CONSTANTS_EP'),
                line,
            )

        if in_use_stmt:
            line = replace_routine_names(line, const_renames)
            if not line.rstrip().endswith('&'):
                in_use_stmt = False

        result.append(line)

    return '\n'.join(result)


def migrate_free_form(source: str, rename_map: dict[str, str],
                      kind: int) -> str:
    """Migrate a free-form Fortran file (.f90).

    These files typically use parameterized types via:
        integer, parameter :: wp = kind(1.d0)  -- or --
        integer, parameter :: wp = real64       (iso_fortran_env)
    We change the wp/sp/dp definition and rename routines.
    """
    source = rewrite_la_constants_use(source, kind)
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        nl = '\n' if line.endswith('\n') else ''

        # Change working-precision parameter definition
        stripped = _replace_kind_parameter(stripped, kind)
        if not stripped.lstrip().startswith('!'):
            stripped = _strip_iso_fortran_env_realN(stripped)
            if not stripped:
                # The USE line was entirely dropped (realN was its only
                # import); skip it so we don't leave a blank gap.
                continue
            stripped = replace_intrinsic_calls(stripped, kind)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, kind)
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


def migrate_file_to_string(src_path: Path, rename_map: dict[str, str],
                           kind: int,
                           parser: str | None = None,
                           parser_cmd: str | None = None,
                           ) -> tuple[str, str] | None:
    """Migrate a file in memory. Returns (target_filename, migrated_text).

    Separated from :func:`migrate_file` so the pipeline can perform a
    convergence check across co-family members that map to the same
    output name before committing to disk.

    Args:
        parser: Which parse tree backend to use: ``'flang'``,
            ``'gfortran'``, or ``None`` (regex-only, no compiler).
        parser_cmd: Explicit path to the compiler binary.  When *None*
            the compiler is located via PATH.
    """
    ext = src_path.suffix.lower()
    source = src_path.read_text(errors='replace')

    facts = None
    if parser == 'flang':
        from .flang_parser import scan_file as flang_scan, find_flang
        cmd = parser_cmd or find_flang()
        if cmd:
            facts = flang_scan(src_path, cmd)
    elif parser == 'gfortran':
        from .gfortran_parser import scan_file as gfortran_scan, find_gfortran
        cmd = parser_cmd or find_gfortran()
        if cmd:
            facts = gfortran_scan(src_path, cmd)

    if facts is not None:
        migrated = _migrate_with_flang(source, ext, rename_map, kind, facts)
    elif ext in ('.f', '.for'):
        migrated = migrate_fixed_form(source, rename_map, kind)
    elif ext in ('.f90', '.f95', '.F90'):
        migrated = migrate_free_form(source, rename_map, kind)
    else:
        return None

    out_name = target_filename(src_path.name, rename_map)
    return out_name, migrated


def migrate_file(src_path: Path, output_dir: Path,
                 rename_map: dict[str, str], kind: int,
                 parser: str | None = None,
                 parser_cmd: str | None = None) -> str | None:
    """Migrate a single Fortran source file. Returns output filename.

    Uses the specified *parser* backend (``'flang'``, ``'gfortran'``,
    or ``None`` for regex-only) to guide transformations.
    """
    result = migrate_file_to_string(src_path, rename_map, kind, parser,
                                    parser_cmd)
    if result is None:
        return None
    out_name, migrated = result
    (output_dir / out_name).write_text(migrated)
    return out_name


def _migrate_with_flang(source: str, ext: str,
                        rename_map: dict[str, str], kind: int,
                        facts) -> str:
    """Migrate using Flang parse tree facts to guide transformations.

    The parse tree tells us WHAT to transform (which types, which names,
    which literals exist in the file). We then apply source-level
    replacements using the same functions as the regex path, but with
    confidence from the parse tree.

    Key benefits over pure regex:
      - Type declarations are identified structurally (no false positives
        from type keywords in comments or string literals on code lines)
      - Intrinsic calls vs. type keywords are disambiguated
      - Call sites are identified precisely
    """
    # Build the effective rename map for this file — only include names
    # that actually appear in the parse tree as routine defs, call sites,
    # or external declarations.
    file_names: set[str] = set()
    for rd in facts.routine_defs:
        file_names.add(rd.name)
    for cs in facts.call_sites:
        file_names.add(cs.name)
    for en in facts.external_names:
        file_names.add(en)

    # Restrict rename map to names present in this file
    file_rename_map = {
        k: v for k, v in rename_map.items()
        if k in file_names
    }

    # Determine which floating-point types exist (from parse tree)
    has_float_types = any(
        td.type_spec in ('DoublePrecision', 'Real', 'Complex')
        for td in facts.type_decls
    )

    has_real_literals = bool(facts.real_literals)

    # Now apply source-level transformations line by line
    if ext in ('.f90', '.f95', '.F90'):
        return _migrate_free_form_flang(
            source, file_rename_map, kind, has_float_types)
    else:
        return _migrate_fixed_form_flang(
            source, file_rename_map, kind,
            has_float_types, has_real_literals)


def _migrate_fixed_form_flang(source: str, rename_map: dict[str, str],
                              kind: int, has_float_types: bool,
                              has_real_literals: bool) -> str:
    """Fixed-form migration guided by Flang parse tree."""
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
            if has_float_types:
                stripped = replace_type_decls(stripped, kind)
        else:
            if has_float_types:
                stripped = replace_type_decls(stripped, kind)
                stripped = replace_standalone_real_complex(stripped, kind)
            if has_real_literals:
                stripped = replace_literals(stripped, kind)
            stripped = replace_intrinsic_calls(stripped, kind)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, kind)
            stripped = replace_routine_names(stripped, rename_map)
            stripped = replace_xerbla_strings(stripped, rename_map)
            stripped = reformat_fixed_line(stripped)
        result.append(stripped + nl)
    return ''.join(result)


def _migrate_free_form_flang(source: str, rename_map: dict[str, str],
                             kind: int, has_float_types: bool) -> str:
    """Free-form migration guided by Flang parse tree."""
    source = rewrite_la_constants_use(source, kind)
    lines = source.splitlines(keepends=True)
    result = []
    for line in lines:
        stripped = line.rstrip('\n\r')
        nl = '\n' if line.endswith('\n') else ''
        stripped = _replace_kind_parameter(stripped, kind)
        if not stripped.lstrip().startswith('!'):
            stripped = _strip_iso_fortran_env_realN(stripped)
            if not stripped:
                continue
            stripped = replace_intrinsic_calls(stripped, kind)
            stripped = replace_intrinsic_decls(stripped)
            stripped = replace_generic_conversions(stripped, kind)
        stripped = replace_routine_names(stripped, rename_map)
        if stripped.lstrip().startswith('!') and has_float_types:
            stripped = replace_type_decls(stripped, kind)
        result.append(stripped + nl)
    return ''.join(result)
