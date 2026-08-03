"""Tests for the two paren-call rewriters in ``intrinsics_rw``.

``replace_intrinsic_calls`` and ``replace_generic_conversions`` share one
scan/splice loop (``_rewrite_paren_calls``) and differ only in the
decider each hands it. The point of these tests is to pin those
differences, because they are the reason the two functions cannot simply
be merged:

* the ``add_kind`` ``KIND=`` skip test is depth-tracking in one and
  depth-blind in the other, so ``REAL(f(x, KIND=1))`` is rewritten by
  one and skipped by the other;
* only the generic rewriter refuses a call that already has its full
  argument count;
* only the intrinsic rewriter has the ambiguous-type-spec guard, and
  only it falls back to a bare rename when the closing paren is off the
  line.

Collapsing either KIND= predicate into the other changes migrated
output; these are the cases that would flag it.
"""

from migrator.fortran.intrinsics_rw import (
    replace_generic_conversions,
    replace_intrinsic_calls,
)
from migrator.target_mode import load_target

KIND16 = load_target('kind16')
MULTIFLOATS = load_target('multifloats')
K = KIND16.kind_suffix


# --- the divergent KIND= predicates ---------------------------------

def test_nested_kind_blocks_the_generic_rewriter():
    """Depth-blind: a ``KIND=`` anywhere inside the parens skips."""
    line = '      X = REAL(F(Y, KIND=1))'
    assert replace_generic_conversions(line, KIND16) == line


def test_nested_kind_does_not_block_the_intrinsic_rewriter():
    """Depth-tracking: only a *top-level* ``KIND=`` skips, so the same
    argument text gets the opposite verdict here."""
    out = replace_intrinsic_calls('      X = DBLE(F(Y, KIND=1))', KIND16)
    assert out == f'      X = REAL(F(Y, KIND=1), KIND={K})'


def test_top_level_kind_blocks_both():
    line = '      X = REAL(Y, KIND=8)'
    assert replace_generic_conversions(line, KIND16) == line
    assert replace_intrinsic_calls('      X = DBLE(Y, KIND=8)', KIND16) == \
        '      X = DBLE(Y, KIND=8)'


# --- guards only one of them has ------------------------------------

def test_generic_rewriter_skips_a_full_argument_list():
    """``REAL`` takes one argument, ``CMPLX`` two; a call already at its
    limit has the kind positionally or is already migrated."""
    assert replace_generic_conversions('      X = REAL(Y, WP)', KIND16) == \
        '      X = REAL(Y, WP)'
    assert replace_generic_conversions('      X = CMPLX(A, B, WP)',
                                       KIND16) == '      X = CMPLX(A, B, WP)'


def test_generic_rewriter_counts_a_trailing_comma():
    """The comma count is a count of commas, not of arguments — a
    trailing one still means the list is full."""
    line = '      X = REAL(Y,)'
    assert replace_generic_conversions(line, KIND16) == line


def test_intrinsic_rewriter_has_no_argument_count_guard():
    """Pinning current behavior, not endorsing it: the intrinsic
    rewriter appends a fourth argument here where the generic one backs
    off (previous test). Reachable only for a 3-argument CMPLX, which
    the corpus does not contain."""
    assert replace_intrinsic_calls('      X = CMPLX(A, B, WP)', KIND16) == \
        f'      X = CMPLX(A, B, WP, KIND={K})'


def test_ambiguous_type_spec_guard_is_intrinsic_only():
    """``CMPLX(WP)`` might be the type spec ``CMPLX(KIND=WP)``, so the
    intrinsic rewriter leaves it alone. The generic rewriter has no such
    guard — its lookbehind already restricts it to expression context,
    where a type specification cannot appear."""
    assert replace_intrinsic_calls('      X = CMPLX(WP)', KIND16) == \
        '      X = CMPLX(WP)'
    assert replace_generic_conversions('      X = CMPLX(WP)', KIND16) == \
        f'      X = CMPLX(WP, KIND={K})'


def test_bare_integer_is_ambiguous_for_the_intrinsic_rewriter_only():
    assert replace_intrinsic_calls('      X = CMPLX(3)', KIND16) == \
        '      X = CMPLX(3)'
    assert replace_generic_conversions('      X = CMPLX(3)', KIND16) == \
        f'      X = CMPLX(3, KIND={K})'


def test_unbalanced_paren_falls_back_to_a_bare_rename():
    """Only the intrinsic rewriter has an ``on_unbalanced`` handler; the
    generic one stops scanning the line."""
    assert replace_intrinsic_calls('      X = DBLE(Y +', KIND16) == \
        '      X = REAL(Y +'
    line = '      X = REAL(Y +'
    assert replace_generic_conversions(line, KIND16) == line


# --- shared loop mechanics ------------------------------------------

def test_every_call_on_the_line_is_rewritten():
    out = replace_intrinsic_calls('      X = DBLE(A) + DBLE(B)', KIND16)
    assert out == f'      X = REAL(A, KIND={K}) + REAL(B, KIND={K})'


def test_scanning_resumes_past_the_replacement():
    """The generic rewriter does not rename, so its output still matches
    its own pattern: resuming anywhere inside the replacement would
    rewrite the rewrite (or spin forever). Only the outermost call is
    touched."""
    out = replace_generic_conversions('      X = REAL(REAL(A))', KIND16)
    assert out == f'      X = REAL(REAL(A), KIND={K})'


def test_continuation_marker_is_restored():
    """``replace_generic_conversions`` masks column 6 while scanning."""
    out = replace_generic_conversions('     +   X = REAL(Y)', KIND16)
    assert out == f'     +   X = REAL(Y, KIND={K})'


# --- the shared multifloats replacement builder ---------------------

def test_multifloats_real_drops_the_kind_argument_in_both():
    """Both routes now go through ``_wrap_constructor_call``; the kind
    spec is dropped because the constructor takes one component."""
    assert replace_intrinsic_calls('      X = DBLE(Y, WP)', MULTIFLOATS) == \
        f'      X = {MULTIFLOATS.real_constructor}(Y)'
    assert replace_generic_conversions('      X = REAL(Y)', MULTIFLOATS) == \
        f'      X = {MULTIFLOATS.real_constructor}(Y)'


def test_multifloats_cmplx_prewraps_its_arguments_in_both():
    rc, cc = MULTIFLOATS.real_constructor, MULTIFLOATS.complex_constructor
    assert replace_intrinsic_calls('      X = DCMPLX(A, B)', MULTIFLOATS) == \
        f'      X = {cc}({rc}(A),{rc}(B))'
    assert replace_generic_conversions('      X = CMPLX(A, B)',
                                       MULTIFLOATS) == \
        f'      X = {cc}({rc}(A),{rc}(B))'
