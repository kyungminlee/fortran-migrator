"""Narrow multifloats (real64x2 / cmplx64x2) items in formatted output
statements so each presents a single value against an ordinary real edit
descriptor.

Background
----------
The migrator retypes a genuine-MUMPS ``DOUBLE PRECISION`` statistic such
as ``RINFOG(40)`` to ``TYPE(real64x2)`` (two ``limbs``). A formatted
``WRITE`` whose format supplies a single real edit descriptor for that
item then expands the derived type into its two components (F2018
13.7.6): the second limb has no descriptor left, triggers format
reversion, and — when reversion lands on an integer descriptor — aborts
at runtime ("Expected INTEGER for item ... got REAL"). Other sites
silently misprint (extra number, shifted columns, spurious records).

Fix: inside formatted output lists only, wrap a *direct reference* to a
real64x2 variable as ``dble(x)`` and a cmplx64x2 reference as
``cmplx(dble(x%re), dble(x%im), kind=8)`` — one intrinsic-real value per
descriptor, matching how the genuine scalar-real arithmetics print. Both
``dble`` and ``cmplx`` are multifloats generics that extend the
intrinsics, so the emitted calls dispatch correctly and, because this
pass runs before USE specialization, ``_build_use_only_clause`` imports
whichever it sees referenced.

Continuation across cpp directives
----------------------------------
The crashing WRITE (``dana_aux.F``) has ``#if``/``#endif`` directives
interleaved in its continuation lines. The fixed-form segmenter breaks a
statement at every preprocessor line, so the fragment that actually
carries ``RINFOG(1)`` — ``     &  KEEP(56), KEEP(61), RINFOG(1)`` — is a
headless continuation with no ``WRITE`` keyword. A single-statement pass
keyed on finding ``WRITE``/``PRINT`` never sees it. The caller therefore
drives this module *statefully*: :func:`narrow_multifloats_io_open`
reports whether a statement opened a formatted output list, and while it
stays open every subsequent fixed-form continuation fragment
(:func:`is_fixed_io_continuation`) is narrowed via
:func:`narrow_io_continuation`. The list is considered closed as soon as
a non-continuation statement begins.

Scope guards
------------
- Only fires when the per-file real/complex name oracle is non-empty —
  i.e. module (multifloats) targets. Kind-based targets (kind10/kind16
  are intrinsic reals) pass empty sets and this is a no-op.
- Only *pure data references* (an identifier with optional subscript,
  chained by ``%``) are narrowed. Expressions, intrinsic calls, literals
  and implied-do lists are left untouched — narrowing an arbitrary
  expression would need value-type inference and risks wrapping an
  integer item in ``dble`` (which would then feed a REAL to an integer
  descriptor — the very crash we fix).
- Unformatted writes (MUMPS save/restore of ``id%RINFO``) are never
  touched — narrowing a binary record would corrupt the saved state.
"""
import re
from dataclasses import dataclass

from .lex import (
    FIXED_FORM_MARGIN,
    _find_inline_bang,
    is_continuation_line,
    match_paren,
    string_mask,
)

_KW_RE = re.compile(r'(write|print)\b', re.IGNORECASE)
_IDENT_RE = re.compile(r'[A-Za-z_]\w*')


@dataclass(frozen=True)
class NameOracle:
    """The four name sets the narrowing passes consult, as one value.

    They only ever travel together — deleting any one makes the rest
    meaningless — and the two operations *on* them (:meth:`is_empty`,
    :meth:`uppercased`) used to be free functions taking all four.
    Bundling them also removes the parameter-rebinding hazard the loose
    form had, where a caller reassigned two of its own four parameters to
    their uppercased forms halfway through the body.

    ``real`` / ``cplx`` are the per-file variable-name oracles;
    ``real_comp`` / ``cplx_comp`` are the global derived-type *component*
    oracles, matched only against ``%``-qualified references.
    """
    real: frozenset[str] = frozenset()
    cplx: frozenset[str] = frozenset()
    real_comp: frozenset[str] = frozenset()
    cplx_comp: frozenset[str] = frozenset()

    def is_empty(self) -> bool:
        """True when every oracle is empty — the no-op signal for
        kind-based targets (kind10/kind16 are intrinsic reals)."""
        return not (self.real or self.cplx or self.real_comp or self.cplx_comp)

    def uppercased(self) -> 'NameOracle':
        """This oracle with every name upper-cased, for case-insensitive
        matching against Fortran source."""
        return NameOracle(
            frozenset(n.upper() for n in self.real),
            frozenset(n.upper() for n in self.cplx),
            frozenset(n.upper() for n in self.real_comp),
            frozenset(n.upper() for n in self.cplx_comp),
        )


def _split_top_level(line: str, inside: list[bool], start: int, end: int):
    """(s, e) spans of comma-separated items in ``line[start:end]``,
    splitting only on commas outside strings and at paren depth 0."""
    spans = []
    s = start
    depth = 0
    for i in range(start, end):
        if inside[i]:
            continue
        c = line[i]
        if c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
        elif c == ',' and depth == 0:
            spans.append((s, i))
            s = i + 1
    spans.append((s, end))
    return spans


def _reference_last_name(text: str) -> tuple[str | None, bool]:
    """If ``text`` is a pure data reference — an identifier with optional
    subscript, chained by ``%`` — return ``(NAME, qualified)`` where NAME is
    the UPPERCASE final component and ``qualified`` is True when the
    reference reached it through at least one ``%`` (i.e. a derived-type
    component access such as ``id%DKEEP(160)``). Otherwise ``(None, False)``
    (expression, literal, intrinsic call, implied-do, ...).
    """
    t = text.strip()
    if not t:
        return None, False
    i, n = 0, len(t)
    last = None
    qualified = False
    while True:
        m = _IDENT_RE.match(t, i)
        if not m:
            return None, False
        last = m.group(0)
        i = m.end()
        if i < n and t[i] == '(':
            # Quote-blind: a subscript expression in a data reference never
            # carries a string literal, and this text has no mask to hand.
            j = match_paren(t, i)
            if j is None:
                return None, False  # unbalanced
            i = j + 1
        while i < n and t[i] == ' ':
            i += 1
        if i >= n:
            return last.upper(), qualified
        if t[i] == '%':
            qualified = True
            i += 1
            while i < n and t[i] == ' ':
                i += 1
            continue
        return None, False  # trailing operator / other -> not a pure reference


def _is_formatted_control(line: str, inside: list[bool],
                          open_idx: int, close_idx: int) -> bool:
    """A WRITE control list ``(...)`` is formatted iff it carries a
    format: an explicit ``FMT=`` keyword, or a positional second item.
    Unformatted writes (unit only, or unit + keyword specs like
    ``IOSTAT=``/``NML=``) return False and are left untouched.
    """
    spans = _split_top_level(line, inside, open_idx + 1, close_idx)
    for (s, e) in spans:
        if re.match(r'fmt\s*=', line[s:e].strip(), re.IGNORECASE):
            return True
    if len(spans) >= 2:
        second = line[spans[1][0]:spans[1][1]].strip()
        # A keyword spec is ``NAME = value`` (single ``=``, not ``==``).
        # Anything else in the second position is the positional format.
        if second and not re.match(r'[A-Za-z]\w*\s*=(?!=)', second):
            return True
    return False


def _wrap(itemtext: str, kind: str) -> str:
    """Wrap the item core, preserving surrounding whitespace for layout."""
    lead = len(itemtext) - len(itemtext.lstrip())
    trail = len(itemtext) - len(itemtext.rstrip())
    core = itemtext.strip()
    if kind == 'real':
        w = f'dble({core})'
    else:
        w = f'cmplx(dble({core}%re), dble({core}%im), kind=8)'
    head = itemtext[:lead]
    tail = itemtext[len(itemtext) - trail:] if trail else ''
    return head + w + tail


def _narrow_items(line: str, inside: list[bool], list_start: int,
                  stmt_end: int, oracle: NameOracle) -> str:
    """Wrap every pure real64x2/cmplx64x2 reference in the output list
    ``line[list_start:stmt_end]``.

    ``oracle`` must already be upper-cased (see
    :meth:`NameOracle.uppercased`). Its ``real`` / ``cplx`` sets are the
    per-file oracle (local variables and dummies); they match a reference
    by its final name regardless of form. Its ``real_comp`` / ``cplx_comp``
    sets are the *global* derived-type component oracle; they match only a
    ``%``-qualified reference (``id%DKEEP(160)``), never a bare name — a
    bare identifier that merely shares a struct field's name is a
    different, locally-typed variable.

    Returns ``line`` unchanged when nothing matched (so the caller's
    ``s == joined`` verbatim path still fires)."""
    if list_start >= stmt_end:
        return line
    spans = _split_top_level(line, inside, list_start, stmt_end)
    out = [line[:list_start]]
    cursor = list_start
    changed = False
    for (s, e) in spans:
        out.append(line[cursor:s])
        item = line[s:e]
        name, qualified = _reference_last_name(item)
        is_real = name is not None and (
            name in oracle.real or (qualified and name in oracle.real_comp))
        is_complex = name is not None and (
            name in oracle.cplx or (qualified and name in oracle.cplx_comp))
        if is_real:
            out.append(_wrap(item, 'real'))
            changed = True
        elif is_complex:
            out.append(_wrap(item, 'complex'))
            changed = True
        else:
            out.append(item)
        cursor = e
    out.append(line[cursor:])
    return ''.join(out) if changed else line


def narrow_multifloats_io_open(line: str, real_names: set[str],
                               complex_names: set[str],
                               real_comp: set[str] = frozenset(),
                               complex_comp: set[str] = frozenset(),
                               ) -> tuple[str, bool]:
    """Narrow real64x2/cmplx64x2 references in a formatted WRITE/PRINT
    output list. Returns ``(new_line, opened)`` where ``opened`` is True
    when ``line`` is a formatted output statement — the signal the caller
    uses to keep narrowing subsequent continuation fragments of the same
    logical statement (see module docstring). No-op (returns
    ``(line, False)``) unless some oracle is non-empty (module targets) and
    the statement is a formatted output statement.

    ``real_comp`` / ``complex_comp`` are the global derived-type component
    oracle — matched only against ``%``-qualified references (see
    :func:`_narrow_items`).
    """
    oracle = NameOracle(frozenset(real_names), frozenset(complex_names),
                        frozenset(real_comp), frozenset(complex_comp))
    if oracle.is_empty():
        return line, False
    if not line:
        return line, False
    low = line.lower()
    if 'write' not in low and 'print' not in low:
        return line, False

    oracle = oracle.uppercased()
    inside = string_mask(line)

    # Statement start: skip an optional numeric label + whitespace.
    pos = re.match(r'\s*(?:\d+\s+)?', line).end()

    # Optional logical-IF prefix: ``IF (cond) <output-stmt>``.
    ifm = re.match(r'if\s*\(', line[pos:], re.IGNORECASE)
    if ifm:
        open_idx = pos + ifm.end() - 1  # the '(' of the condition
        if inside[open_idx]:
            return line, False
        close_idx = match_paren(line, open_idx, inside=inside)
        if close_idx is None:
            return line, False
        pos = close_idx + 1
        while pos < len(line) and line[pos] == ' ':
            pos += 1

    km = _KW_RE.match(line, pos)
    if not km:
        return line, False
    kw = km.group(1).lower()
    after = km.end()

    stmt_end = _find_inline_bang(line)  # len(line) when no inline comment

    if kw == 'write':
        pm = re.match(r'\s*\(', line[after:])
        if not pm:
            return line, False  # not a WRITE(...) — e.g. ``write = value``
        open_idx = after + pm.end() - 1
        if inside[open_idx]:
            return line, False
        close_idx = match_paren(line, open_idx, inside=inside)
        if close_idx is None:
            return line, False
        if not _is_formatted_control(line, inside, open_idx, close_idx):
            return line, False  # unformatted write — never narrow
        list_start = close_idx + 1
    else:  # print — format is mandatory, so always formatted
        list_start = None
        depth = 0
        for i in range(after, stmt_end):
            if inside[i]:
                continue
            c = line[i]
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            elif c == ',' and depth == 0:
                list_start = i + 1
                break
        if list_start is None:
            return line, False  # PRINT with no output list

    new = _narrow_items(line, inside, list_start, stmt_end, oracle)
    return new, True


def narrow_multifloats_io(line: str, real_names: set[str],
                          complex_names: set[str],
                          real_comp: set[str] = frozenset(),
                          complex_comp: set[str] = frozenset()) -> str:
    """Single-statement narrowing — thin wrapper over
    :func:`narrow_multifloats_io_open` that drops the ``opened`` flag.
    Used by the free-form path and by callers that never need to track a
    cpp-interrupted continuation."""
    return narrow_multifloats_io_open(
        line, real_names, complex_names, real_comp, complex_comp)[0]


def is_fixed_io_continuation(line: str) -> bool:
    """True when ``line`` is a fixed-form continuation line (a non-blank,
    non-``0`` character in column 6, columns 1–5 blank). Mirrors the
    continuation test in ``_segment_fixed_form_statements`` so a
    cpp-split fragment is recognised the same way the segmenter split it
    (in particular a tab in column 6 is NOT a continuation marker).
    """
    return is_continuation_line(line, tab_marker=False)


def narrow_io_continuation(line: str, real_names: set[str],
                           complex_names: set[str],
                           real_comp: set[str] = frozenset(),
                           complex_comp: set[str] = frozenset()) -> str:
    """Narrow every pure real64x2/cmplx64x2 reference in a fixed-form
    continuation fragment whose parent formatted output list is still
    open. The column-6 continuation marker (index 5) is preserved; the
    output-list items begin at column 7 (index 6)."""
    oracle = NameOracle(frozenset(real_names), frozenset(complex_names),
                        frozenset(real_comp), frozenset(complex_comp))
    if oracle.is_empty():
        return line
    if not is_fixed_io_continuation(line):
        return line
    inside = string_mask(line)
    stmt_end = _find_inline_bang(line)
    return _narrow_items(line, inside, FIXED_FORM_MARGIN, stmt_end,
                         oracle.uppercased())
