"""Unit tests for the multifloats branch of pyengine.fortran_migrator.

These tests cover the per-pass transformations in isolation plus a few
small end-to-end migrations of synthetic BLAS-style sources. They run
without any external compiler — they only check the textual output of
the migrator. The expensive end-to-end smoke test that compiles real
BLAS sources against the multifloats module lives outside pytest.
"""

import textwrap

import pytest

from pyengine.target_mode import multifloats_target, kind_target
from pyengine.fortran_migrator import (
    replace_type_decls,
    replace_literals,
    replace_known_constants,
    strip_known_constants_from_decls,
    convert_parameter_stmts,
    convert_data_stmts,
    insert_use_multifloats,
    rewrite_la_constants_use,
    _la_constants_rename_map,
    migrate_fixed_form,
    migrate_free_form,
)


# ---------------------------------------------------------------------------
# TargetMode factory
# ---------------------------------------------------------------------------


def test_multifloats_target_basic_shape():
    mf = multifloats_target()
    assert mf.name == 'multifloats'
    assert mf.is_kind_based is False
    assert mf.real_type == 'TYPE(float64x2)'
    assert mf.complex_type == 'TYPE(complex128x2)'
    assert mf.real_constructor == 'float64x2'
    assert mf.complex_constructor == 'complex128x2'
    assert mf.module_name == 'multifloats'
    # W (real) / U (complex) prefix per MULTIFLOATS.md
    assert mf.prefix_map['R'] == 'W'
    assert mf.prefix_map['C'] == 'U'


def test_kind_target_is_kind_based():
    k16 = kind_target(16)
    assert k16.is_kind_based is True
    assert k16.real_type == 'REAL(KIND=16)'
    assert k16.kind_suffix == 16


def test_multifloats_known_constants_excludes_rtmin_rtmax():
    """RTMIN/RTMAX collide with local LAPACK variables — must NOT be in
    the auto-rename set."""
    mf = multifloats_target()
    assert 'RTMIN' not in mf.known_constants
    assert 'RTMAX' not in mf.known_constants
    assert 'ZERO' in mf.known_constants
    assert 'ONE' in mf.known_constants


# ---------------------------------------------------------------------------
# replace_type_decls
# ---------------------------------------------------------------------------


@pytest.fixture
def mf():
    return multifloats_target()


@pytest.mark.parametrize('src,expected', [
    ('      DOUBLE PRECISION X', '      TYPE(float64x2) X'),
    ('      DOUBLE COMPLEX Z', '      TYPE(complex128x2) Z'),
    ('      COMPLEX*16 Z', '      TYPE(complex128x2) Z'),
    ('      REAL*8 X', '      TYPE(float64x2) X'),
    ('      REAL*4 X', '      TYPE(float64x2) X'),
])
def test_replace_type_decls_keywords(mf, src, expected):
    assert replace_type_decls(src, mf) == expected


def test_replace_type_decls_preserves_attrs(mf):
    src = '      DOUBLE PRECISION, INTENT(IN) :: A(LDA,*)'
    out = replace_type_decls(src, mf)
    assert 'TYPE(float64x2)' in out
    assert 'INTENT(IN)' in out
    assert 'A(LDA,*)' in out


# ---------------------------------------------------------------------------
# strip_known_constants_from_decls
# ---------------------------------------------------------------------------


def test_strip_drops_entire_line_when_all_names_known(mf):
    src = '      DOUBLE PRECISION ZERO,ONE\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert new == ''
    assert removed == {'ZERO': 'MF_ZERO', 'ONE': 'MF_ONE'}


def test_strip_filters_partial_decl(mf):
    src = '      DOUBLE PRECISION DFLAG,TWO,W,Z,ZERO\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert 'DFLAG' in new
    assert 'W' in new and 'Z' in new
    assert 'TWO' not in new
    assert 'ZERO' not in new
    assert removed == {'TWO': 'MF_TWO', 'ZERO': 'MF_ZERO'}


def test_strip_handles_multiline_continuation(mf):
    """Continuation line carries names that should also be filtered."""
    src = (
        '      DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP,\n'
        '     $                 DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO\n'
    )
    new, removed = strip_known_constants_from_decls(src, mf)
    assert 'ONE' not in new
    assert 'TWO' not in new
    assert 'ZERO' not in new
    assert 'GAM' in new
    assert 'GAMSQ' in new
    assert 'RGAMSQ' in new
    assert removed == {'ONE': 'MF_ONE', 'TWO': 'MF_TWO', 'ZERO': 'MF_ZERO'}


def test_strip_matches_bare_real_keyword(mf):
    """sgemm.f uses bare REAL (single precision) — must be matched too."""
    src = '      REAL ONE,ZERO\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert new == ''
    assert removed == {'ONE': 'MF_ONE', 'ZERO': 'MF_ZERO'}


def test_strip_does_not_touch_array_specs(mf):
    """A line with array specs is not safe to filter — leave it alone."""
    src = '      DOUBLE PRECISION A(LDA,*),X(*),Y(*)\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert new == src
    assert removed == {}


def test_strip_skips_function_return_decl(mf):
    src = '      DOUBLE PRECISION FUNCTION FOO(N)\n'
    new, removed = strip_known_constants_from_decls(src, mf)
    assert new == src
    assert removed == {}


def test_strip_kind_target_is_noop():
    src = '      DOUBLE PRECISION ZERO,ONE\n'
    new, removed = strip_known_constants_from_decls(src, kind_target(16))
    assert new == src
    assert removed == {}


# ---------------------------------------------------------------------------
# replace_literals
# ---------------------------------------------------------------------------


def test_replace_literals_d_exponent(mf):
    src = "      X = 1.0D+0 + 2.5D-3"
    out = replace_literals(src, mf)
    assert "float64x2('1.0D+0')" in out
    assert "float64x2('2.5D-3')" in out


def test_replace_literals_bare_float(mf):
    src = "      X = 0.5"
    out = replace_literals(src, mf)
    assert "float64x2('0.5D0')" in out


def test_replace_literals_wp_suffix(mf):
    src = "   x = 0.5_wp + 1.0_wp"
    out = replace_literals(src, mf)
    assert "float64x2('0.5D0')" in out
    assert "float64x2('1.0D0')" in out


def test_replace_literals_does_not_double_wrap(mf):
    """Idempotence: re-running on already-wrapped output is a no-op."""
    src = "      X = float64x2('1.0D+0')"
    out1 = replace_literals(src, mf)
    out2 = replace_literals(out1, mf)
    # Second pass mustn't add another constructor layer.
    assert out1.count('float64x2(') == 1
    assert out2.count('float64x2(') == 1


def test_replace_literals_preserves_strings(mf):
    """Float-shaped text inside string literals must be untouched."""
    src = "      WRITE(*,*) '1.0D+0 in a string'"
    out = replace_literals(src, mf)
    assert "'1.0D+0 in a string'" in out
    assert 'float64x2(' not in out


def test_replace_literals_skips_fortran_operators(mf):
    """``.EQ.`` etc. must not be parsed as floating-point literals."""
    src = "      IF (X.EQ.1.0D+0) GO TO 10"
    out = replace_literals(src, mf)
    assert '.EQ.' in out
    assert "float64x2('1.0D+0')" in out


# ---------------------------------------------------------------------------
# replace_known_constants
# ---------------------------------------------------------------------------


def test_replace_known_constants_uses_per_file_renames(mf):
    src = "      X = ZERO + ONE"
    out = replace_known_constants(src, mf, renames={'ZERO': 'MF_ZERO', 'ONE': 'MF_ONE'})
    assert out == "      X = MF_ZERO + MF_ONE"


def test_replace_known_constants_skips_comment_lines(mf):
    src = "*     performs ONE pass through ZERO array"
    out = replace_known_constants(src, mf, renames={'ZERO': 'MF_ZERO', 'ONE': 'MF_ONE'})
    assert out == src


def test_replace_known_constants_skips_inline_comment_text(mf):
    src = "      X = ZERO     ! the ZERO constant"
    out = replace_known_constants(src, mf, renames={'ZERO': 'MF_ZERO'})
    # ZERO before the ! gets renamed; ZERO after stays as English.
    assert out == "      X = MF_ZERO     ! the ZERO constant"


def test_replace_known_constants_skips_string_literals(mf):
    src = "      CALL XERBLA('ZERO INPUT')"
    out = replace_known_constants(src, mf, renames={'ZERO': 'MF_ZERO'})
    assert out == src


def test_replace_known_constants_skips_decl_lines(mf):
    src = "      INTEGER ZERO"
    out = replace_known_constants(src, mf, renames={'ZERO': 'MF_ZERO'})
    assert out == src


def test_replace_known_constants_kind_mode_is_noop():
    src = "      X = ZERO"
    out = replace_known_constants(src, kind_target(16), renames={'ZERO': 'MF_ZERO'})
    assert out == src


# ---------------------------------------------------------------------------
# convert_parameter_stmts / convert_data_stmts
# ---------------------------------------------------------------------------


def test_convert_parameter_drops_known_constant(mf):
    src = '      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)\n'
    new, assigns, dropped = convert_parameter_stmts(src, mf)
    # Both names are known constants — no assignments emitted, line elided.
    assert assigns == []
    assert dropped == {'ONE': 'MF_ONE', 'ZERO': 'MF_ZERO'}
    assert 'PARAMETER' not in new or 'Converted' in new


def test_convert_parameter_keeps_unknown_fp(mf):
    src = '      PARAMETER (R=ONE/IPW2)\n'
    new, assigns, dropped = convert_parameter_stmts(src, mf)
    assert any('R = ONE/IPW2' in a for a in assigns)
    assert dropped == {}


def test_convert_parameter_mixed_type_keeps_integers(mf):
    src = '      PARAMETER (LV=128,IPW2=4096,R=1.0D0)\n'
    new, assigns, dropped = convert_parameter_stmts(src, mf)
    # Integers stay in the rewritten PARAMETER stmt, FP becomes assignment
    assert 'LV=128' in new
    assert 'IPW2=4096' in new
    assert 'PARAMETER' in new
    assert any('R = 1.0D0' in a for a in assigns)


def test_convert_data_drops_known_constants(mf):
    src = '      DATA ZERO,TWO/0.D0,2.D0/\n'
    new, assigns, dropped = convert_data_stmts(src, mf)
    assert assigns == []
    assert dropped == {'ZERO': 'MF_ZERO', 'TWO': 'MF_TWO'}


def test_convert_data_keeps_unknown_names(mf):
    src = '      DATA GAM,GAMSQ/4096.D0,16777216.D0/\n'
    new, assigns, dropped = convert_data_stmts(src, mf)
    assert any('GAM = 4096.D0' in a for a in assigns)
    assert any('GAMSQ = 16777216.D0' in a for a in assigns)
    assert dropped == {}


# ---------------------------------------------------------------------------
# insert_use_multifloats
# ---------------------------------------------------------------------------


def test_insert_use_multifloats_placement_after_decls(mf):
    """Extra param/data assignments must come AFTER the declaration block,
    not right after the SUBROUTINE header."""
    src = textwrap.dedent('''\
              SUBROUTINE FOO(A)
              DOUBLE PRECISION A
              INTEGER N
              A = 1.0D0
              END
        ''')
    out = insert_use_multifloats(src, mf, extra_lines=['      X = float64x2(\'1.0D0\')\n'])
    lines = out.splitlines()
    use_idx = next(i for i, l in enumerate(lines) if 'USE multifloats' in l)
    assign_idx = next(i for i, l in enumerate(lines) if 'X = float64x2' in l)
    int_idx = next(i for i, l in enumerate(lines) if 'INTEGER N' in l)
    real_idx = next(i for i, l in enumerate(lines) if 'DOUBLE PRECISION A' in l)
    # USE comes right after the header, before any declaration
    assert use_idx < real_idx
    # Extra assignment comes AFTER the declaration block
    assert assign_idx > int_idx


def test_insert_use_multifloats_dedupes(mf):
    src = textwrap.dedent('''\
              SUBROUTINE FOO()
              USE multifloats
              END
        ''')
    out = insert_use_multifloats(src, mf, extra_lines=None)
    assert out.count('USE multifloats') == 1


# ---------------------------------------------------------------------------
# _la_constants_rename_map / rewrite_la_constants_use
# ---------------------------------------------------------------------------


def test_la_constants_rename_map_uses_w_u_for_multifloats(mf):
    m = _la_constants_rename_map(mf)
    assert m['DZERO'] == 'WZERO'
    assert m['DSAFMIN'] == 'WSAFMIN'
    assert m['ZZERO'] == 'UZERO'
    # Unprefixed forms must NOT be in the map (would clobber USE alias LHS)
    assert 'ZERO' not in m
    assert 'SAFMIN' not in m


def test_la_constants_rename_map_kind16():
    m = _la_constants_rename_map(kind_target(16))
    assert m['DZERO'] == 'QZERO'
    assert m['ZZERO'] == 'XZERO'


def test_rewrite_la_constants_use_pattern_b(mf):
    src = textwrap.dedent('''\
        subroutine foo
           use LA_CONSTANTS, &
              only: wp=>dp, zero=>dzero, half=>dhalf, one=>done, &
                    safmin=>dsafmin, safmax=>dsafmax
           real(wp) :: x
        end
    ''')
    out = rewrite_la_constants_use(src, mf)
    assert 'LA_CONSTANTS_MF' in out
    # wp=>dp removed, but local aliases preserved with W-prefixed RHS
    assert 'wp=>dp' not in out
    assert 'zero=>wzero' in out
    assert 'safmin=>wsafmin' in out
    assert 'safmax=>wsafmax' in out


# ---------------------------------------------------------------------------
# End-to-end fixed-form smoke
# ---------------------------------------------------------------------------


SYNTHETIC_BLAS = """\
      SUBROUTINE FOO(N,X,ALPHA)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
      IF (ALPHA.EQ.ZERO) RETURN
      DO 10 I = 1,N
          X(I) = ALPHA*X(I) + ONE
   10 CONTINUE
      RETURN
      END
"""


def test_end_to_end_fixed_form(mf):
    out = migrate_fixed_form(SYNTHETIC_BLAS, {}, mf)

    # Header is unchanged
    assert 'SUBROUTINE FOO' in out
    # USE multifloats inserted right after the header
    assert 'USE multifloats' in out
    # Type decls converted
    assert 'TYPE(float64x2) ALPHA' in out
    assert 'TYPE(float64x2) X(*)' in out
    # PARAMETER block + decl for ONE/ZERO are gone (the PARAMETER line
    # may survive as a `! Converted to assignments: ...` comment).
    assert 'DOUBLE PRECISION ONE,ZERO' not in out
    code_lines = [
        l for l in out.splitlines()
        if l and not l.lstrip().startswith('!') and not l.startswith(('*', 'C', 'c'))
    ]
    code = '\n'.join(code_lines)
    assert 'PARAMETER (ONE=' not in code
    # Body references renamed
    assert 'IF (ALPHA.EQ.MF_ZERO) RETURN' in out
    assert '+ MF_ONE' in out
    # No bare ZERO/ONE references in body code
    body_lines = [
        l for l in out.splitlines()
        if l and not l.startswith(('*', 'C', 'c', '!')) and 'Converted' not in l
    ]
    body = '\n'.join(body_lines)
    import re as _re
    assert not _re.search(r'\bZERO\b', body)
    assert not _re.search(r'\bONE\b', body)


SYNTHETIC_LAPACK_FREE_FORM = """\
subroutine FOO( x, y )
   use LA_CONSTANTS, &
      only: wp=>dp, zero=>dzero, one=>done, safmin=>dsafmin
   real(wp) :: x, y
   if( x == zero ) then
      y = one / safmin
   end if
end subroutine
"""


def test_end_to_end_free_form_pattern_b(mf):
    out = migrate_free_form(SYNTHETIC_LAPACK_FREE_FORM, {}, mf)
    assert 'LA_CONSTANTS_MF' in out
    assert 'wp=>dp' not in out
    # Local aliases preserved (lowercase) with W-prefixed RHS
    assert 'zero=>wzero' in out
    assert 'safmin=>wsafmin' in out
    # real(wp) → TYPE(float64x2)
    assert 'TYPE(float64x2)' in out
    # Body references stay as the local alias names
    assert 'x == zero' in out
    assert 'one / safmin' in out
