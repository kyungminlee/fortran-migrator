"""PARAMETER / DATA statement conversion (Cluster H).

Retypes constants declared in PARAMETER and DATA statements to the target
kind. Extracted verbatim from ``fortran_migrator.py``.
"""
import re

from ..target_mode import TargetMode
from .lex import (
    _is_fp_value, _join_continued_lines, _scope_indices,
    split_top_level_commas,
)
from .decls import _scan_complex_var_names


def _is_complex_value(name: str, val: str, complex_names: set[str],
                      target_mode: TargetMode, *,
                      type_is_complex: bool = False) -> bool:
    """True when a PARAMETER entry carries complex semantics: the value
    is a complex literal / constructor call, the local declaration of
    the name is COMPLEX, or (combined form) the declared type-spec is
    complex. Such names must NOT be folded into the multifloats
    real-constant rename map — they stay runtime assignments so the
    variable retains its complex type."""
    cx_ctor = (target_mode.complex_constructor or '').lower()
    return (
        type_is_complex
        or ('(' in val and ',' in val)
        or bool(cx_ctor and cx_ctor in val.lower())
        or 'cmplx' in val.lower()
        or 'dcmplx' in val.lower()
        or name.upper() in complex_names
    )


class CppGuardStack:
    """The cpp conditionals in effect at the current source line.

    Feed it every line in order; ask it to :meth:`wrap` any runtime
    assignment it emits. The assignment then stays in scope exactly when
    the declaration it replaces was in scope — otherwise gfortran sees an
    assignment to an undeclared symbol.

    Each frame remembers two things: the ``#if`` line to re-emit, and the
    history of conditions already seen in this if-block. The history is
    what makes ``#elif`` / ``#else`` work — an assignment in the ``#else``
    branch must be wrapped in the *negation* of the branches above it,
    not in the (unrelated) opening ``#if``.
    """

    def __init__(self) -> None:
        self._frames: list[dict] = []

    @staticmethod
    def _condition(directive: str) -> str:
        s = directive.strip()
        if s.startswith('#ifdef'):
            return f'defined({s[len("#ifdef"):].strip()})'
        if s.startswith('#ifndef'):
            return f'!defined({s[len("#ifndef"):].strip()})'
        if s.startswith('#elif'):
            return s[len('#elif'):].strip()
        if s.startswith('#if'):
            return s[len('#if'):].strip()
        return ''

    def feed(self, line: str) -> None:
        """Update the stack for one source line (a no-op unless it is a
        conditional directive)."""
        s = line.lstrip()
        if s.startswith('#if'):
            self._frames.append({
                'wrap': line.rstrip('\n'),
                'history': [self._condition(line)],
            })
        elif s.startswith('#elif'):
            if self._frames:
                frame = self._frames[-1]
                negated = ' || '.join(f'({c})' for c in frame['history'])
                new_cond = self._condition(line)
                frame['history'].append(new_cond)
                frame['wrap'] = f'#if !({negated}) && ({new_cond})'
        elif s.startswith('#else'):
            if self._frames:
                frame = self._frames[-1]
                negated = ' || '.join(f'({c})' for c in frame['history'])
                frame['wrap'] = f'#if !({negated})'
        elif s.startswith('#endif'):
            if self._frames:
                self._frames.pop()

    def wrap(self, assn: str) -> str:
        """Wrap ``assn`` in the effective ``#if``/``#endif`` guards."""
        if not self._frames:
            return assn
        pre = ''.join(f['wrap'] + '\n' for f in self._frames)
        return pre + assn + '#endif\n' * len(self._frames)


def _split_param_entry(entry: str) -> tuple[str, str] | None:
    """``'ONE = 1.0D+0'`` -> ``('ONE', '1.0D+0')``; ``None`` when the
    entry carries no initializer."""
    if '=' not in entry:
        return None
    name, val = entry.split('=', 1)
    return name.strip(), val.strip()


def _has_complex_literal_entry(entries: list[str]) -> bool:
    """True when some entry's value looks like a complex literal —
    a parenthesised, comma-separated pair."""
    return any(
        ('(' in v.split('=', 1)[1] and ',' in v.split('=', 1)[1])
        for v in entries if '=' in v
    )


# Combined declaration-attribute form: ``TYPE, PARAMETER :: name = val``.
# The type-spec captured in group 2 is whatever comes before the
# ``, PARAMETER ::`` token. After matching, the type-spec is preserved
# for the new declaration line and each ``name = val`` entry is split off
# as a runtime assignment — same logic as the standalone-PARAMETER branch.
_COMBINED_PARAM_RE = re.compile(
    r'^(\s+)('
    r'(?:DOUBLE\s+PRECISION|DOUBLE\s+COMPLEX'
    r'|REAL\s*\*\s*\d+|COMPLEX\s*\*\s*\d+'
    r'|REAL\s*\(\s*[^)]*\)|COMPLEX\s*\(\s*[^)]*\)'
    r'|TYPE\s*\(\s*[A-Za-z_]\w*\s*\)'
    r'|INTEGER\s*\(\s*[^)]*\)|INTEGER'
    r'|REAL|COMPLEX)'
    r')\s*,\s*PARAMETER\s*::\s*(.+?)\s*(!.*)?$',
    re.IGNORECASE,
)

_PARAM_RE = re.compile(
    r'^(\s{6,}|^\s*)PARAMETER\s*\((.*)\)\s*(!.*)?$', re.IGNORECASE)

# Fixed-form statement opener, used only to decide whether joining
# continuation lines is worth trying.
_PARAM_STMT_START_RE = re.compile(r'^\s{6,}PARAMETER\b', re.IGNORECASE)

_INTEGER_TYPE_SPEC_RE = re.compile(r'\s*INTEGER', re.IGNORECASE)
_COMPLEX_TYPE_SPEC_RE = re.compile(r'COMPLEX|cmplx64x2', re.IGNORECASE)


def _convert_combined_param(
    cm: re.Match, target_mode: TargetMode, complex_names: set[str],
    guards: CppGuardStack,
) -> tuple[str, list[str], dict[str, str]] | None:
    """Convert one ``TYPE, PARAMETER :: name = val [, ...]`` statement.

    Splits it into (a) a plain declaration line listing only the variable
    names and (b) a runtime assignment per entry. Known-constant names
    with non-complex values are dropped entirely — the multifloats module
    supplies the constant.

    Returns ``(decl_line, assignments, dropped_known)``, where
    ``decl_line`` is ``''`` when every entry was dropped, or ``None``
    when the statement must be left exactly as it stands.
    """
    indent, type_sp, items, comment = (
        cm.group(1), cm.group(2), cm.group(3), cm.group(4) or '')
    entries = split_top_level_commas(items, strip=True)
    type_is_complex = bool(_COMPLEX_TYPE_SPEC_RE.search(type_sp))

    # Standard INTEGER PARAMETER lines (``INTEGER, PARAMETER :: N = 5``
    # / ``wp = kind(1.d0)`` etc.) are constant expressions and must stay
    # as PARAMETERs — converting them to runtime assignments would break
    # things like the ``kind(1.d0)`` idiom (the literal gets wrapped to a
    # derived-type value the ``kind()`` intrinsic rejects).
    #
    # The exception is the upstream-bug case
    # ``INTEGER, PARAMETER :: ZERO = (0.0D0, 0.0D0)``
    # (zsol_fwd_aux.F:1095): gfortran tolerates these via numeric
    # coercion, but multifloats mode rewrites the literal to
    # ``cmplx64x2(...)`` which is not a constant expression for INTEGER.
    # When we detect a complex-literal value on an INTEGER LHS, override
    # the declared type to the multifloats complex type so the runtime
    # assignment is type-correct (the use sites are downstream complex
    # contexts in every observed case).
    if _INTEGER_TYPE_SPEC_RE.match(type_sp):
        if not _has_complex_literal_entry(entries):
            return None
        type_sp = target_mode.complex_type
        type_is_complex = True

    kept_names: list[str] = []
    assignments: list[str] = []
    dropped: dict[str, str] = {}
    for entry in entries:
        split = _split_param_entry(entry)
        if split is None:
            # A PARAMETER without a value is not legal Fortran, so this
            # branch should not trigger in practice.
            return None
        name, val = split
        if not _is_fp_value(val, target_mode.known_constants):
            # INTEGER PARAMETER, character literal, etc. don't need
            # conversion — leave the combined-form line alone.
            return None

        is_cx_value = _is_complex_value(
            name, val, complex_names, target_mode,
            type_is_complex=type_is_complex,
        )
        if name.upper() in target_mode.known_constants and not is_cx_value:
            dropped[name.upper()] = target_mode.known_constants[name.upper()]
            continue

        kept_names.append(name)
        assignments.append(guards.wrap(f"{indent}{name} = {val}{comment}\n"))

    decl = (f"{indent}{type_sp} :: {', '.join(kept_names)}{comment}\n"
            if kept_names else '')
    return decl, assignments, dropped


def _convert_standalone_param(
    m: re.Match, joined: str, target_mode: TargetMode,
    complex_names: set[str], guards: CppGuardStack,
) -> tuple[str, list[str], dict[str, str]]:
    """Convert one standalone ``PARAMETER (name = val, ...)`` statement.

    Returns ``(replacement_line, assignments, dropped_known)``. The
    replacement is ``''`` when every entry was a known constant supplied
    by the multifloats module, so the line simply disappears.
    """
    indent, params_content, comment = m.group(1), m.group(2), m.group(3) or ''

    kept_parts: list[str] = []
    assignments: list[str] = []
    dropped: dict[str, str] = {}
    for part in split_top_level_commas(params_content):
        split = _split_param_entry(part)
        if split is None:
            kept_parts.append(part)
            continue
        name, val = split
        if not _is_fp_value(val, target_mode.known_constants):
            kept_parts.append(part)
            continue
        is_cx_value = _is_complex_value(name, val, complex_names, target_mode)
        if name.upper() in target_mode.known_constants and not is_cx_value:
            dropped[name.upper()] = target_mode.known_constants[name.upper()]
            continue
        assignments.append(guards.wrap(f"{indent}{name} = {val}{comment}\n"))

    if kept_parts:
        replacement = f"{indent}PARAMETER ({', '.join(kept_parts)}){comment}\n"
    elif assignments:
        # Some FP entries became runtime assignments — leave a short
        # marker comment so reviewers can find the source.
        replacement = f"{indent}! Converted to assignments below: {joined.strip()}\n"
    else:
        replacement = ''
    return replacement, assignments, dropped


def convert_parameter_stmts(
    source: str, target_mode: TargetMode,
) -> tuple[str, list[tuple[int, str]], dict[str, str]]:
    """Convert floating-point PARAMETER statements to executable assignments.

    Returns ``(new_source, fp_assignments, dropped_known)`` where
    ``fp_assignments`` is a list of ``(scope_index, assignment_text)``
    tuples so the caller can insert each assignment into the correct
    procedure scope. ``dropped_known`` maps each known-constant name
    skipped from the PARAMETER list to its multifloats replacement (so
    the caller can add it to the per-file rename set).

    Multi-line PARAMETER statements (fixed-form column-6 continuation)
    are joined into a single logical statement before parsing. The
    original line(s) are replaced as a unit so the line count of the
    output may differ from the input.
    """
    if target_mode.is_kind_based:
        return source, [], {}

    # Pre-scan declarations so we can tell whether a name like ``ONE``
    # was declared COMPLEX. The original LAPACK convention in
    # Z-prefixed routines is ``COMPLEX*16 ONE; PARAMETER(ONE = 1.0D+0)``
    # — the value is a real literal but it carries complex semantics
    # because Fortran promotes the literal to the declared type. The
    # multifloats migrator must NOT fold such names into the real
    # ``MF_ONE`` constant.
    complex_names = _scan_complex_var_names(source)

    lines = source.splitlines(keepends=True)
    scope_vec = _scope_indices(lines)
    result, fp_assignments = [], []
    dropped_known: dict[str, str] = {}
    guards = CppGuardStack()

    i = 0
    while i < len(lines):
        line = lines[i]
        guards.feed(line)

        # Try matching a single-line PARAMETER first; if not, try
        # joining continuation lines and matching the joined form.
        joined, next_i = line.rstrip('\n'), i + 1
        if (_PARAM_RE.match(joined) is None
                and _PARAM_STMT_START_RE.match(joined)):
            joined, next_i = _join_continued_lines(lines, i)

        cm = _COMBINED_PARAM_RE.match(joined)
        if cm:
            converted = _convert_combined_param(
                cm, target_mode, complex_names, guards)
            if converted is not None:
                decl, assignments, dropped = converted
                fp_assignments.extend((scope_vec[i], a) for a in assignments)
                dropped_known.update(dropped)
                if decl:
                    result.append(decl)
                i = next_i
                continue
            # Bailed out — emit the line unchanged. The standalone branch
            # below cannot match a combined form, so fall straight
            # through to the unchanged-emit path.
            result.append(line)
            i += 1
            continue

        m = _PARAM_RE.match(joined)
        if m:
            replacement, assignments, dropped = _convert_standalone_param(
                m, joined, target_mode, complex_names, guards)
            fp_assignments.extend((scope_vec[i], a) for a in assignments)
            dropped_known.update(dropped)
            if replacement:
                result.append(replacement)
            i = next_i
            continue

        result.append(line)
        i += 1
    return "".join(result), fp_assignments, dropped_known


def convert_data_stmts(
    source: str, target_mode: TargetMode,
) -> tuple[str, list[tuple[int, str]], dict[str, str]]:
    """Convert floating-point DATA statements to executable assignments.

    Returns ``(new_source, fp_assignments, dropped_known)`` — see
    :func:`convert_parameter_stmts` for the meaning of the third tuple
    element.
    """
    if target_mode.is_kind_based:
        return source, [], {}

    lines = source.splitlines(keepends=True)
    scope_vec = _scope_indices(lines)
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
            tmp_vals = [v.strip() for v in split_top_level_commas(vals_part)]
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
                    scope = scope_vec[i]
                    fp_assignments.extend((scope, a) for a in line_assignments)
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


def _param_nuke_map(target_mode: TargetMode) -> dict[str, str]:
    """Names of free-form file-scope PARAMETER declarations the target
    module supplies itself, mapped to the module's spelling.

    ``la_constants_map`` is the same list of constants one layer up (the
    ``USE LA_CONSTANTS_MW`` rename map); the free-form Pattern A files
    ``USE multifloats`` directly but want the same names dropped, so it
    is the base here too. ``param_nuke_extra`` adds the names that are
    not LA_CONSTANTS members at all — see the target YAML.

    Both are lower-cased: the caller matches against a lower-cased line.
    """
    return {k.lower(): v
            for k, v in (*target_mode.la_constants_map.items(),
                         *target_mode.param_nuke_extra.items())}


def nuke_multifloats_params(source: str, removed_known: dict[str, str],
                            target_mode: TargetMode) -> str:
    """Comment out free-form PARAMETER declarations the target supplies.

    Heuristic pass over free-form source: declaration lines that mention
    a name from :func:`_param_nuke_map` are commented out (including
    ``&``-continued follow-on lines), and each dropped name is recorded
    in *removed_known* (mutated in place) so
    ``replace_known_constants`` later rewrites its uses to the ``DD_*``
    constant.
    """
    nuke_renames = _param_nuke_map(target_mode)
    lines_tmp = source.splitlines()
    res_tmp = []
    in_comment_block = False
    nuke_names = set(nuke_renames)

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
    return '\n'.join(res_tmp)
