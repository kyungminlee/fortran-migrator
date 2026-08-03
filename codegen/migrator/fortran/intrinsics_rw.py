"""Intrinsic call / declaration / conversion rewriting (Cluster C).

Rewrites intrinsic calls and generic type conversions to their migrated
equivalents via the INTRINSIC_* maps. Extracted verbatim from
``fortran_migrator.py``.
"""
import re

from ..intrinsics import INTRINSIC_MAP, INTRINSIC_DECL_MAP
from ..target_mode import TargetMode
from .lex import (
    FIXED_FORM_LABEL_FIELD, FIXED_FORM_MARGIN, FIXED_FORM_WIDTH,
    is_continuation_line, is_maskable_continuation, match_paren,
    split_top_level_commas,
)

# The INTRINSIC declaration, recognised once. ``replace_intrinsic_decls``
# and ``_dedup_intrinsic_stmts`` both gate on the first and both capture
# the keyword prefix with the second. (``strip_overloaded_intrinsics``
# uses a third, multiline pattern of its own — it runs on a joined
# logical line and accepts a zero-width indent, so it is deliberately not
# one of these.)
_INTRINSIC_STMT_RE = re.compile(r'\s+INTRINSIC\b', re.IGNORECASE)
_INTRINSIC_PREFIX_RE = re.compile(r'(\s+INTRINSIC\s+)', re.IGNORECASE)

# Continuation prefix for a re-wrapped INTRINSIC statement: column-6 '$'
# marker, then the padding that lines the names up under the first
# chunk's. Distinct from ``fixedform._CONT_PREFIX`` ('+' and no padding)
# — the migrated corpus records both spellings byte for byte.
_INTRINSIC_CONT_PREFIX = FIXED_FORM_LABEL_FIELD + '$' + ' ' * 19


# ``REAL(`` / ``CMPLX(`` paren-call detector for ``replace_generic_conversions``.
# The name is group 1: the pattern also consumes the whitespace before it,
# so ``m.start()`` is not where the name begins (see ``name_group``).
_GENERIC_CONV_RE: dict[str, re.Pattern] = {
    _name: re.compile(rf'(?<=[=+\-*/,(.\0])\s*\b({_name})\s*\(', re.IGNORECASE)
    for _name in ('REAL', 'CMPLX')
}


# Precompiled once per process — INTRINSIC_MAP keys are static so the
# per-name patterns never vary. ``replace_intrinsic_calls`` is on the
# hot path (called per statement from migrate_*_form), so rebuilding
# 74 patterns per call was wasteful.
_INTRINSIC_CALL_RE: dict[str, re.Pattern] = {
    _name: re.compile(rf'\b{_name}\s*\(', re.IGNORECASE)
    for _name in INTRINSIC_MAP
}


_INTRINSIC_CALL_RE_REPL: dict[str, re.Pattern] = {
    _name: re.compile(rf'\b({_name})(\s*\()', re.IGNORECASE)
    for _name in INTRINSIC_MAP
}


# Cheap early-out gate: if no intrinsic name from the map appears
# followed by an opening paren on this line, the 74-iteration body
# below has nothing to do. One alternation search per line is much
# cheaper than 74 single-name searches.
_INTRINSIC_CALL_GATE_RE = re.compile(
    r'\b(?:' + '|'.join(INTRINSIC_MAP) + r')\s*\(',
    re.IGNORECASE,
)


def _is_ambiguous_type_spec(old_name: str, inner: str) -> bool:
    """True when ``old_name(inner)`` cannot be told apart from a type
    specification and must therefore be left alone.

    Bare integer literals and bare symbols are skipped only for the
    ambiguous names — ``REAL(3)`` might be the type spec ``REAL(KIND=3)``.
    ``DBLE(3)`` / ``DCMPLX(3)`` / etc. are unambiguously conversion-function
    calls and must be rewritten (otherwise an ``s*`` sibling that wrote
    ``real(3)`` diverges post-migration: ``real(3)`` strips to ``3`` in the
    light-diff normalizer; ``dble(3)`` does not).
    """
    inner_stripped = inner.strip().upper()
    is_type_spec_name = old_name.upper() in ('REAL', 'CMPLX', 'COMPLEX')
    return bool(
        re.match(r'KIND\s*=', inner_stripped)
        or (is_type_spec_name and inner_stripped.isdigit())
        or (is_type_spec_name and re.match(r'^[A-Z_]\w*$', inner_stripped))
    )


def _has_top_level_kind(inner: str) -> bool:
    """True when ``inner`` already carries a ``KIND=`` argument at paren
    depth 0 — adding a second one would be a duplicate keyword argument."""
    depth = 0
    for ii, ch in enumerate(inner):
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        elif (depth == 0 and inner[ii:ii + 4].upper() == 'KIND'
                and re.match(r'KIND\s*=', inner[ii:], re.IGNORECASE)):
            return True
    return False


def _wrap_constructor_call(
    old_name: str, new_name: str, inner: str, target_mode: TargetMode,
    real_names: set[str] | None, complex_names: set[str] | None,
) -> str:
    """Build the ``wrap_constructor`` (multifloats) replacement text for
    ``old_name(inner)``.

    ``DBLE`` / ``DREAL`` / ``REAL`` map to the ``float64x2(...)`` generic
    constructor, which handles every input type uniformly: integer,
    real(sp/dp), float64x2 (identity), and complex (extract real part).
    Older code routed these through ``MF_REAL`` but that interface only
    accepts float64x2 / complex128x2, breaking ``DBLE(integer)`` from
    LAPACK.

    ``CMPLX`` / ``DCMPLX`` map to ``complex128x2``. multifloats's
    complex128x2 interface only has overloads for (float64x2[, float64x2])
    and (real(dp)[, real(dp)]) — there's no integer overload, so a 1-arg
    call with an integer argument falls back to the structure constructor
    and fails. :func:`_wrap_complex_args` pre-wraps any single-arg cmplx
    with float64x2(...) so it always picks the cx_from_mf_1 procedure.

    ``AIMAG`` / ``CONJG`` are module generics — leave the name alone.
    """
    if old_name.upper() in ('REAL', 'DREAL', 'DBLE'):
        # Drop the optional kind-spec second argument when present:
        # `REAL(x, wp)` / `REAL(x, KIND=wp)` / `DBLE(x, wp)` becomes
        # `real64x2(x)` for multifloats. The constructor takes a single
        # value argument; passing the kind-symbol as a second positional
        # component overflows the 1-component derived type ("Too many
        # components in structure constructor"). Split at the top-level
        # comma respecting nested parens.
        args = split_top_level_commas(inner)
        inner_first = args[0] if args else inner
        return f'{target_mode.real_constructor}({inner_first})'
    if old_name.upper() in ('CMPLX', 'DCMPLX'):
        wrapped = _wrap_complex_args(
            inner, target_mode, real_names, complex_names=complex_names,
        )
        return f'{target_mode.complex_constructor}({wrapped})'
    return f'{new_name}({inner})'


def _rewrite_paren_calls(
    line: str, pattern: re.Pattern, decide, *, name_group: int = 0,
    on_unbalanced=None,
) -> str:
    """Rewrite every ``name(...)`` call ``pattern`` finds on ``line``.

    The scan/splice mechanics shared by :func:`replace_intrinsic_calls`
    and :func:`replace_generic_conversions`: find the name, match its
    closing paren, hand the argument text to ``decide``, splice the
    answer in and resume *past the replacement* — never inside it, or a
    rewritten call would be rescanned and rewritten again.

    ``decide(inner)`` returns the replacement text for the whole
    ``name(inner)`` span, or ``None`` to leave this call alone and skip
    to the next. Everything the two callers disagree about — which skip
    heuristics apply, which ``KIND=`` test, how the multifloats
    constructor is spelled — lives in their deciders, not here.

    ``name_group`` is the capture group holding the name, i.e. where the
    replacement starts. ``_INTRINSIC_CALL_RE`` has no groups and begins
    at the name (0); ``_GENERIC_CONV_RE`` swallows the leading
    whitespace and captures the name in group 1.

    ``on_unbalanced(line, m)`` handles a call whose closing paren is not
    on this line, returning ``(line, search_start)``. The default is to
    stop scanning: no closing paren means no argument text, and every
    later match on the line would hit the same wall.
    """
    search_start = 0
    while True:
        m = pattern.search(line, search_start)
        if not m:
            return line
        start = m.start(name_group)
        paren_start = line.index('(', start)
        close_pos = match_paren(line, paren_start)

        if close_pos is None:
            if on_unbalanced is None:
                return line
            line, search_start = on_unbalanced(line, m)
            continue

        replacement = decide(line[paren_start + 1:close_pos])
        if replacement is None:
            search_start = close_pos + 1
            continue

        line = line[:start] + replacement + line[close_pos + 1:]
        search_start = start + len(replacement)


def _intrinsic_call_decider(
    old_name: str, new_name: str, target_mode: TargetMode,
    real_names: set[str] | None, complex_names: set[str] | None,
):
    """Decider for one ``INTRINSIC_MAP`` entry whose ``needs_kind`` flag
    is set — a call whose argument list has to be inspected, not just
    renamed.

    Unlike :func:`_generic_conversion_decider` this has no
    argument-count guard, and deliberately so: that decider only adds a
    kind, so backing off a full argument list is free, while this one
    also *renames*. Skipping ``DBLE(a, b, wp)`` would leave ``DBLE`` in
    a kind16 build — silently double precision, where appending the
    kind at least fails to compile. Handling it properly means a third
    behavior (drop the positional kind, as the ``wrap_constructor`` arm
    already does), and no such call exists in the corpus, so neither is
    implemented. ``test_intrinsic_rewriter_has_no_argument_count_guard``
    pins the current answer as a tripwire.
    """
    def decide(inner: str) -> str | None:
        if _is_ambiguous_type_spec(old_name, inner):
            return None
        if target_mode.intrinsic_mode == 'add_kind':
            if _has_top_level_kind(inner):
                return None
            return f'{new_name}({inner}, KIND={target_mode.kind_suffix})'
        return _wrap_constructor_call(
            old_name, new_name, inner, target_mode, real_names,
            complex_names,
        )
    return decide


def _rename_bare_intrinsic(old_name: str, new_name: str):
    """``on_unbalanced`` handler for :func:`replace_intrinsic_calls`.

    The call is split across a continuation we cannot see, so there is
    no argument list to inspect: fall back to renaming the bare name and
    keep scanning the rest of the line.
    """
    def handle(line: str, m: re.Match) -> tuple[str, int]:
        if old_name.upper() == new_name.upper():
            return line, m.end()
        start = m.start()
        matched = line[start:m.end() - 1]
        repl = new_name.upper() if matched.isupper() else new_name.lower()
        return (line[:start] + repl + line[start + len(matched):],
                start + len(repl))
    return handle


def replace_intrinsic_calls(
    line: str,
    target_mode: TargetMode,
    real_names: set[str] | None = None,
    complex_names: set[str] | None = None,
) -> str:
    """Replace type-specific intrinsic function calls."""
    if not _INTRINSIC_CALL_GATE_RE.search(line):
        return line
    for old_name, (new_name, needs_kind) in INTRINSIC_MAP.items():
        if needs_kind:
            line = _rewrite_paren_calls(
                line, _INTRINSIC_CALL_RE[old_name],
                _intrinsic_call_decider(old_name, new_name, target_mode,
                                        real_names, complex_names),
                on_unbalanced=_rename_bare_intrinsic(old_name, new_name),
            )
        else:
            def _call_replace(m, _new=new_name):
                matched_name = m.group(1)
                rest = m.group(2)
                return (_new.upper() if matched_name.isupper() else _new.lower()) + rest

            line = _INTRINSIC_CALL_RE_REPL[old_name].sub(_call_replace, line)
    return line


def _generic_conversion_decider(name: str, target_mode: TargetMode,
                                complex_names: set[str] | None):
    """Decider for a generic ``REAL(`` / ``CMPLX(`` call in expression
    context.

    What makes this *not* a spelling of :func:`_intrinsic_call_decider`:
    a call already carrying its full argument count is left alone
    (``REAL`` takes 1, ``CMPLX`` 2), since the kind argument would be a
    third and ``REAL(x, KIND=wp)`` is already migrated. The intrinsic
    rewriter deliberately has no such guard — see
    :func:`_intrinsic_call_decider`.

    There is no ``_is_ambiguous_type_spec`` guard either: the pattern's
    lookbehind already restricts the match to expression context, where
    a type specification cannot appear.

    The ``KIND=`` skip test *was* a second, depth-blind predicate
    (a plain ``\\bKIND\\s*=`` search) until it was unified onto
    :func:`_has_top_level_kind`. The two disagreed on
    ``REAL(f(x, KIND=1))``, where the nested ``KIND=`` belongs to ``f``
    and not to the ``REAL`` being rewritten: the depth-blind test
    skipped a call that should be migrated. Staged output is unchanged
    across all five targets — no such call exists in the corpus.
    """
    max_args = 1 if name == 'REAL' else 2

    def decide(inner: str) -> str | None:
        if target_mode.intrinsic_mode != 'add_kind':
            return _wrap_constructor_call(name, name, inner, target_mode,
                                          None, complex_names)
        if _has_top_level_kind(inner):
            return None
        # ``keep_trailing`` makes the piece count a faithful top-level
        # comma count for a trailing comma too ("x," -> 1, not 0).
        if len(split_top_level_commas(inner, keep_trailing=True)) - 1 >= max_args:
            return None
        return f'{name}({inner}, KIND={target_mode.kind_suffix})'
    return decide


def replace_generic_conversions(
    line: str,
    target_mode: TargetMode,
    complex_names: set[str] | None = None,
) -> str:
    """Add KIND (or wrap in constructor) to generic REAL() and CMPLX() calls in expression context."""
    # Deliberately stricter than lex.is_continuation_line: comment
    # markers ('!', 'C', 'c', '*') and tab in col 6 must NOT be masked
    # with '\0' below, or the marker restore would corrupt them.
    is_fixed_cont = is_maskable_continuation(line)
    cont_marker = line[5] if is_fixed_cont else ''
    if is_fixed_cont:
        line = line[:5] + '\0' + line[6:]

    for name in ('REAL', 'CMPLX'):
        line = _rewrite_paren_calls(
            line, _GENERIC_CONV_RE[name],
            _generic_conversion_decider(name, target_mode, complex_names),
            name_group=1,
        )
    if is_fixed_cont:
        line = line[:5] + cont_marker + line[6:]
    return line


def _split_dcolon_sep(text: str) -> tuple[str, str]:
    """Split an optional leading ``::`` separator (with surrounding
    whitespace) off ``text``. Returns ``(sep, rest)`` where ``sep`` is
    the verbatim separator text ('' when absent)."""
    sep_m = re.match(r'\s*::\s*', text)
    if sep_m:
        return text[:sep_m.end()], text[sep_m.end():]
    return '', text


def _normalize_intrinsic_names(names_text: str,
                               target_mode: TargetMode | None) -> list[str]:
    """Comma-split an INTRINSIC name list, dedupe case-insensitively
    (order-preserving, first spelling wins) and, in wrap_constructor
    mode, drop names the target module overloads — gfortran refuses a
    USE-associated generic name when an INTRINSIC declaration in the
    same scope binds it to the standard intrinsic of incompatible
    signature."""
    names = [n.strip() for n in names_text.split(',') if n.strip()]
    seen: set[str] = set()
    deduped: list[str] = []
    for n in names:
        key = n.upper()
        if key not in seen:
            seen.add(key)
            deduped.append(n)
    if target_mode is not None and target_mode.intrinsic_mode == 'wrap_constructor':
        overloaded = frozenset(n.upper() for n in target_mode.module_generic_names)
        deduped = [n for n in deduped if n.upper() not in overloaded]
    return deduped


def replace_intrinsic_decls(line: str, target_mode: TargetMode | None = None) -> str:
    """Replace intrinsic names in INTRINSIC declarations.

    In multifloats mode, additionally drop names that are overloaded by
    the multifloats module — gfortran refuses to accept a USE-associated
    generic name when an INTRINSIC declaration in the same scope binds
    the name to the standard intrinsic of incompatible signature.
    """
    if not _INTRINSIC_STMT_RE.match(line):
        return line
    for old_name, new_name in INTRINSIC_DECL_MAP.items():
        line = re.sub(rf'\b{old_name}\b', new_name, line, flags=re.IGNORECASE)

    m = _INTRINSIC_PREFIX_RE.match(line)
    if m:
        # Name list is this physical line only — a trailing newline is
        # re-attached below, and nothing past it belongs to this statement.
        prefix, name_list = m.group(1), line[m.end():].partition('\n')[0]
        newline = '\n' if line.endswith('\n') else ''
        stripped = name_list.rstrip().rstrip('\n')
        trail = ''
        if stripped.endswith('&'):
            trail = ' &'
            stripped = stripped[:-1]
        elif stripped.endswith(','):
            trail = ','
            stripped = stripped[:-1]
        # Strip an optional ``::`` separator after the keyword.
        sep, stripped = _split_dcolon_sep(stripped)
        deduped = _normalize_intrinsic_names(stripped, target_mode)
        if not deduped:
            # Whole declaration empty — drop the line so we don't emit a
            # bare ``INTRINSIC`` keyword. The trailing comment, if any,
            # is preserved as a regular comment line.
            return ''
        line = prefix + sep + ', '.join(deduped) + trail + newline
    return line


def _collect_intrinsic_stmt(lines: list[str], i: int) -> tuple[list[str], int]:
    """Gather the physical lines of the (possibly continued) INTRINSIC
    statement starting at ``lines[i]``, returning ``(stmt_lines, next_i)``."""
    stmt_lines = [lines[i]]
    j = i + 1
    while j < len(lines):
        next_line = lines[j]
        if is_continuation_line(next_line):
            stmt_lines.append(next_line)
            j += 1
        elif stmt_lines[-1].rstrip().endswith('&'):
            stmt_lines.append(next_line)
            j += 1
        else:
            break
    return stmt_lines, j


def _wrap_intrinsic_stmt(prefix: str, names: list[str]) -> list[str]:
    """Render ``prefix`` + ``names`` as one or more physical lines, none
    wider than the fixed-form limit. Every line but the last ends in a
    comma; every line but the first starts with the continuation prefix."""
    full = prefix + ', '.join(names)
    if len(full) <= FIXED_FORM_WIDTH:
        return [full]
    # -1 leaves room for the trailing comma on every chunk but the last.
    first_cap = FIXED_FORM_WIDTH - len(prefix) - 1
    cont_cap = FIXED_FORM_WIDTH - len(_INTRINSIC_CONT_PREFIX) - 1
    chunks: list[str] = []
    cur = ''
    for name in names:
        addition = (', ' + name) if cur else name
        cap = first_cap if not chunks else cont_cap
        if not cur: cur = name
        elif len(cur) + len(addition) <= cap: cur += addition
        else:
            chunks.append(cur)
            cur = name
    if cur: chunks.append(cur)
    out = []
    for ci, chunk in enumerate(chunks):
        trail = ',' if ci < len(chunks) - 1 else ''
        out.append((prefix if ci == 0 else _INTRINSIC_CONT_PREFIX) + chunk + trail)
    return out


def _dedup_intrinsic_stmts(text: str, target_mode: TargetMode | None = None) -> str:
    """Remove duplicate names from multi-line INTRINSIC statements.

    In multifloats mode, also drop names that are overloaded by the
    multifloats module — they conflict with the use-associated generics.
    """
    lines = text.split('\n')
    result: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if not _INTRINSIC_STMT_RE.match(line):
            result.append(line)
            i += 1
            continue

        stmt_lines, j = _collect_intrinsic_stmt(lines, i)

        m = _INTRINSIC_PREFIX_RE.match(stmt_lines[0])
        if not m:
            result.extend(stmt_lines)
            i = j
            continue

        prefix = m.group(1)
        name_part = stmt_lines[0][len(prefix):]
        for sl in stmt_lines[1:]:
            if is_continuation_line(sl):
                name_part += ' ' + sl[FIXED_FORM_MARGIN:]
            else:
                stripped = sl.lstrip()
                if stripped.startswith('&'): stripped = stripped[1:]
                name_part += ' ' + stripped

        name_part = name_part.replace('&', ' ')
        # Drop an optional ``::`` separator following the keyword.
        _, name_part = _split_dcolon_sep(name_part)
        deduped = _normalize_intrinsic_names(name_part, target_mode)

        if not deduped:
            # Whole INTRINSIC statement empty: drop the entire (possibly
            # multi-line) declaration.
            i = j
            continue

        result.extend(_wrap_intrinsic_stmt(prefix, deduped))
        i = j
    return '\n'.join(result)


def _wrap_complex_args(
    inner: str, target_mode: TargetMode, real_names: set[str] | None,
    complex_names: set[str] | None = None,
) -> str:
    """Pre-wrap each top-level argument of a CMPLX-style call so the
    multifloats complex128x2 interface can dispatch.

    multifloats's ``complex128x2`` interface only has overloads for
    (float64x2[, float64x2]) and (real(dp)[, real(dp)]). For LAPACK's
    common patterns ``CMPLX(N)`` and ``CMPLX(N, 0)`` with integer
    arguments, the call falls back to the structure constructor and
    fails. Wrapping each arg with ``float64x2(...)`` redirects to
    ``cx_from_mf_*``.

    Skip the wrap when the arg is already a float64x2 expression
    (bare ``MF_*`` constant, ``float64x2(...)`` call, or a known
    float64x2 local) — wrapping again would fail because float64x2
    has no identity constructor.

    Also skip the wrap when the arg is itself a complex expression
    (any token in ``complex_names``). ``CMPLX(D, KIND=N)`` /
    ``DCMPLX(D)`` on a complex ``D`` is the Fortran identity (preserves
    Re/Im). Wrapping with ``float64x2(...)`` would dispatch to
    ``dd_from_cdd`` which discards the imaginary part — silently
    corrupting the result.
    """
    parts = split_top_level_commas(inner)

    # Drop a trailing ``KIND=...`` argument. The kind selector is only
    # meaningful in the original CMPLX call signature; the multifloats
    # ``complex128x2`` interface has no equivalent.
    parts = [p for p in parts if not re.match(r'\s*KIND\s*=', p, re.IGNORECASE)]

    # Type-conversion intrinsics that replace_generic_conversions() will
    # handle later — don't pre-wrap them or we get double wrapping.
    _CONVERSION_INTRINSICS = frozenset({
        'REAL', 'DBLE', 'SNGL', 'DREAL', 'DFLOAT', 'FLOAT',
        'CMPLX', 'DCMPLX',
    })

    def _wrap_one(arg: str) -> str:
        s = arg.strip()
        if not s:
            return arg
        # Strip a leading unary +/- and surrounding whitespace.
        body = s.lstrip('+-').lstrip()
        head = re.match(r'([A-Za-z_]\w*)', body)
        # Detect any float64x2 token in the expression — if any operand
        # is a known float64x2 / MF_* / float64x2() call, the whole
        # expression is float64x2 (multifloats provides operator
        # overloads for *, /, +, -). Wrapping the whole expression
        # again with float64x2(...) would fail.
        if body.lower().startswith(target_mode.real_constructor.lower() + '('):
            return arg
        # Skip type-conversion intrinsics — they will be replaced by
        # replace_generic_conversions() later in the pipeline.
        if head and head.group(1).upper() in _CONVERSION_INTRINSICS:
            after_name = body[head.end():]
            if after_name.lstrip().startswith('('):
                return arg
        if real_names:
            for tok in re.finditer(r'\b([A-Za-z_]\w*)\b', body):
                u = tok.group(1).upper()
                if u in real_names or u.startswith('MF_'):
                    return arg
        if complex_names:
            for tok in re.finditer(r'\b([A-Za-z_]\w*)\b', body):
                if tok.group(1).upper() in complex_names:
                    return arg
        if head and head.group(1).upper().startswith('MF_'):
            return arg
        return f'{target_mode.real_constructor}({s})'

    return ','.join(_wrap_one(p) for p in parts)


# Type-conversion intrinsics frequently drift asymmetrically between
# co-family halves: D/Z-half upstream tends to declare ``REAL`` /
# ``CMPLX`` / ``DBLE`` / ``DIMAG`` / ``DCMPLX`` / ``DCONJG`` in its
# INTRINSIC list (used to convert ``DOUBLE COMPLEX`` <-> ``DOUBLE
# PRECISION``), while S/C-half doesn't need them (kind4 default
# handles the same conversions implicitly). After kind16 migration
# both halves use the kind-promoted versions; the asymmetric
# INTRINSIC declarations remain as cosmetic text drift.
# :func:`strip_overloaded_intrinsics` strips this subset on every
# target to converge.
_TYPE_CONV_INTRINSICS = frozenset({
    'REAL', 'CMPLX', 'DCMPLX', 'DBLE',
    'AIMAG', 'CONJG', 'DCONJG', 'DIMAG',
})

# Multifloats also strips the full generic-overload set so INTRINSIC
# declarations of names that ``USE multifloats`` provides as generic
# interfaces don't clash (gfortran: "Cannot change attributes of
# USE-associated symbol").
_MF_OVERLOADED_INTRINSICS = _TYPE_CONV_INTRINSICS | {
    # Unary real -> real
    'ABS', 'SQRT', 'SIN', 'COS', 'TAN', 'EXP', 'LOG', 'LOG10',
    'ATAN', 'ASIN', 'ACOS', 'AINT', 'ANINT', 'SINH', 'COSH',
    'TANH', 'ASINH', 'ACOSH', 'ATANH', 'ERF', 'ERFC',
    'ERFC_SCALED', 'GAMMA', 'LOG_GAMMA', 'BESSEL_J0',
    'BESSEL_J1', 'BESSEL_Y0', 'BESSEL_Y1', 'FRACTION',
    'RRSPACING', 'SPACING', 'EPSILON', 'HUGE', 'TINY',
    # Binary real
    'SIGN', 'MOD', 'ATAN2', 'DIM', 'MODULO', 'HYPOT', 'NEAREST',
    # Variadic
    'MAX', 'MIN',
    # Extra type-conversion
    'INT', 'NINT', 'CEILING', 'FLOOR',
}


def strip_overloaded_intrinsics(line: str, target_mode: TargetMode) -> str:
    """Strip target-supplied names from an INTRINSIC declaration statement.

    Runs on a single joined logical line (via ``_apply_local_passes``),
    so a multi-line INTRINSIC list is cleaned as a whole. Kind-based
    targets strip only the type-conversion subset; multifloats strips
    the full generic-overload set. Returns ``''`` when the kept list is
    empty so the statement is dropped.
    """
    strip_set = (_TYPE_CONV_INTRINSICS if target_mode.is_kind_based
                 else _MF_OVERLOADED_INTRINSICS)

    def clean_intrinsic(m):
        indent, sep, funcs_str, newline = m.group(1), m.group(2), m.group(3), m.group(4)
        kept = [f.strip() for f in funcs_str.split(',')
                if f.strip() and f.strip().upper() not in strip_set]
        return f"{indent}INTRINSIC{sep}{', '.join(kept)}{newline}" if kept else ""

    return re.sub(
        r'(?im)^([ \t]*)INTRINSIC(\s*::\s*|\s+)([A-Za-z0-9_,\s]+?)(\r?\n|$)',
        clean_intrinsic, line,
    )
