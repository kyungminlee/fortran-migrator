"""Fixed-form line reformatting + statement segmentation.

Re-splits migrated fixed-form lines that overflow column 72 and segments
fixed-form statements at their continuation boundaries. Extracted verbatim
from ``fortran_migrator.py``.
"""
from .lex import (
    FIXED_FORM_CODE_WIDTH, FIXED_FORM_LABEL_FIELD, FIXED_FORM_MARGIN,
    FIXED_FORM_WIDTH, _build_split_mask, _ends_in_string, _find_inline_bang,
    is_comment_line, is_continuation_line, split_term,
)

# Column 6 marker emitted on every continuation line this module produces.
CONT_CHAR = '+'
_CONT_PREFIX = FIXED_FORM_LABEL_FIELD + CONT_CHAR

# Split-point search window, in body-relative columns. We prefer to break
# at a comma or space; the search walks back from the last usable column
# but never further than ``_SPLIT_WINDOW`` characters, and never past
# ``_SPLIT_FLOOR`` (which keeps a pathologically token-free line from
# being split down to a stub).
_LAST_CODE_COL = FIXED_FORM_CODE_WIDTH - 1     # == 65
_SPLIT_WINDOW = 30
_SPLIT_FLOOR = 35


def reformat_fixed_line(line: str) -> str:
    # Preprocessor directives (``#if``, ``#include``, ``#define`` ...) are
    # not bound by fixed-form column 72 and must not be split into
    # continuation lines — doing so produces a truncated directive on
    # the first line (e.g. ``#if A || B ||`` dangling) and a second line
    # the preprocessor doesn't understand. Leave them alone regardless
    # of length.
    if line.lstrip().startswith('#'):
        return line
    if (len(line) <= FIXED_FORM_WIDTH or is_comment_line(line)
            or (len(line) > FIXED_FORM_MARGIN
                and line[FIXED_FORM_MARGIN:].lstrip().startswith('!'))):
        return line
    # If an inline ``!`` comment sits within the first 72 columns, keep
    # the whole line intact — fixed-form Fortran ignores columns past
    # 72, and we must NOT split across the comment (the text after
    # ``!`` would otherwise land on a continuation line as code).
    # Scan for the first ``!`` outside a string literal. If it sits
    # within col 72 we must keep the line intact — splitting it would
    # land the comment text on a continuation chunk.
    bang = _find_inline_bang(line)
    if bang < len(line) and bang < FIXED_FORM_WIDTH:
        return line
    prefix, body = (line[:FIXED_FORM_MARGIN] if len(line) >= FIXED_FORM_MARGIN
                    else line.ljust(FIXED_FORM_MARGIN)), line[FIXED_FORM_MARGIN:]
    safe = _build_split_mask(body)
    chunks = []
    while len(body) > FIXED_FORM_CODE_WIDTH:
        split_pos = FIXED_FORM_CODE_WIDTH
        for i in range(_LAST_CODE_COL,
                       max(_SPLIT_FLOOR, _LAST_CODE_COL - _SPLIT_WINDOW), -1):
            if body[i] in (',', ' ') and safe[i]:
                split_pos = i + 1
                break
        else:
            for i in range(_LAST_CODE_COL, 0, -1):
                if safe[i]:
                    split_pos = i
                    break
        chunks.append(body[:split_pos])
        body, safe = body[split_pos:], safe[split_pos:]
    chunks.append(body)
    result_lines = [prefix + chunks[0]]
    for chunk in chunks[1:]: result_lines.append(_CONT_PREFIX + chunk)
    return '\n'.join(result_lines)


def _strip_inline_comment(text: str) -> str:
    """Strip an inline ``!`` comment from a Fortran code line, respecting
    string literals. The trailing whitespace before the ``!`` is also
    trimmed.

    Used by :func:`_segment_fixed_form_statements` when joining
    continuation lines: an inline ``!`` mid-statement would otherwise
    swallow every continuation that follows once
    :func:`reformat_fixed_line` re-splits the joined string at column
    66 and the ``!`` lands on a chunk that includes content from the
    next continuation line. The comment is irretrievably lost from the
    joined / reformatted output, which is the price of correctness.
    Single-physical-line statements with ``s == joined`` go through the
    no-transform path and keep their comments verbatim from ``lines``.
    """
    bang = _find_inline_bang(text)
    if bang < len(text):
        return text[:bang].rstrip()
    return text


def _segment_fixed_form_statements(
    physical: list[str],
) -> list[tuple[str, list[str], list[str], str]]:
    """Group physical fixed-form lines into logical statements.

    Each entry is ``(kind, lines, terminators, joined)`` where ``kind`` is
    ``'blank' | 'comment' | 'pp' | 'code'``. ``lines`` and ``terminators``
    are aligned slices of ``physical`` (text without newline / the
    original line terminator). For ``'code'`` statements with continuation
    lines, ``joined`` is the head plus each continuation's column-7+
    content concatenated with single spaces — with any inline ``!``
    comments stripped (see :func:`_strip_inline_comment`). This is what
    the per-line transform passes operate on, so paren-walkers can match
    across the physical line break. For other kinds, ``joined`` is the
    head text (with inline comments preserved — single-line statements
    don't risk the swallow-the-next-continuation failure mode).
    """
    out: list[tuple[str, list[str], list[str], str]] = []
    i = 0
    while i < len(physical):
        text, term = split_term(physical[i])
        if not text.strip():
            out.append(('blank', [text], [term], text))
            i += 1
            continue
        if text[0] in 'Cc*!':
            out.append(('comment', [text], [term], text))
            i += 1
            continue
        if text.lstrip().startswith('#'):
            out.append(('pp', [text], [term], text))
            i += 1
            continue
        # Code head — absorb any immediately-following continuation lines.
        lines = [text]
        terms = [term]
        joined = text
        j = i + 1
        while j < len(physical):
            ntext, nterm = split_term(physical[j])
            if is_continuation_line(ntext, tab_marker=False):
                lines.append(ntext)
                terms.append(nterm)
                # Strip the previous segment's inline ``!`` comment
                # before appending the next continuation, otherwise the
                # comment swallows every following chunk once
                # ``reformat_fixed_line`` re-splits the joined string.
                # On the first continuation, this also strips any inline
                # comment from the head line.
                head_code = _strip_inline_comment(joined)
                # A character literal split across the continuation must be
                # joined *tight*: in fixed form the string content is the
                # head's text followed directly by the continuation's
                # column-7+ text, with no blank between them. Inserting the
                # token-boundary space here would corrupt the literal (e.g.
                # a WRITE message gaining a doubled space). Outside a string
                # the space is the harmless stand-in for the line break.
                if _ends_in_string(head_code):
                    joined = head_code + ntext[FIXED_FORM_MARGIN:]
                else:
                    joined = head_code + ' ' + ntext[FIXED_FORM_MARGIN:]
                j += 1
            else:
                break
        # Strip any inline comment from the trailing continuation too
        # (so the rebuilt single-line / reformatted output isn't truncated
        # at the comment when reformat_fixed_line walks it).
        if len(lines) > 1:
            joined = _strip_inline_comment(joined)
        out.append(('code', lines, terms, joined))
        i = j
    return out
