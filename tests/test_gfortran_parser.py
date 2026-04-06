"""Tests for pyengine.gfortran_parser.

These tests exercise the ``parse_tree_facts`` function using synthetic
gfortran dump text, so they run without gfortran installed.  An
additional integration test uses real gfortran (skipped if unavailable).
"""

import textwrap
from pathlib import Path

import pytest

from pyengine.gfortran_parser import (
    find_gfortran,
    parse_tree_facts,
    scan_file,
)


# ---------------------------------------------------------------------------
# Synthetic dump for a simple DOUBLE PRECISION subroutine (like DGEMM)
# ---------------------------------------------------------------------------

DGEMM_DUMP = textwrap.dedent("""\
    Namespace: A-H: (REAL 4) I-N: (INTEGER 4) O-Z: (REAL 4)
    procedure name = dgemm
      symtree: 'a'           || symbol: 'a'
        type spec : (REAL 8)
        attributes: (VARIABLE  DIMENSION DUMMY)
        Array spec:(2 [0] AS_ASSUMED_SIZE 1 dgemm:lda 1 () )
      symtree: 'alpha'       || symbol: 'alpha'
        type spec : (REAL 8)
        attributes: (VARIABLE  DUMMY)
      symtree: 'dgemm'       || symbol: 'dgemm'
        type spec : (UNKNOWN 0)
        attributes: (PROCEDURE  SUBROUTINE ARRAY-OUTER-DEPENDENCY)
        Formal arglist: transa transb m n k alpha a lda b ldb beta c ldc
      symtree: 'info'        || symbol: 'info'
        type spec : (INTEGER 4)
        attributes: (VARIABLE )
      symtree: 'lsame'       || symbol: 'lsame'
        type spec : (LOGICAL 4)
        attributes: (PROCEDURE EXTERNAL-PROC  EXTERNAL FUNCTION ARRAY-OUTER-DEPENDENCY)
        result: lsame
      symtree: 'max'         || symbol: 'max'
        type spec : (INTEGER 4)
        attributes: (PROCEDURE INTRINSIC-PROC  INTRINSIC FUNCTION IMPLICIT-TYPE ELEMENTAL PURE)
        result: max
      symtree: 'one'         || symbol: 'one'
        type spec : (REAL 8)
        attributes: (PARAMETER IMPLICIT-SAVE)
        value: 1.0000000000000000e0_8
      symtree: 'temp'        || symbol: 'temp'
        type spec : (REAL 8)
        attributes: (VARIABLE )
      symtree: 'xerbla'      || symbol: 'xerbla'
        type spec : (UNKNOWN 0)
        attributes: (PROCEDURE EXTERNAL-PROC  EXTERNAL SUBROUTINE ARRAY-OUTER-DEPENDENCY)

      code:
      ASSIGN dgemm:nota lsame[[((dgemm:transa) ('N'))]]
      CALL xerbla (('DGEMM ') (dgemm:info))
      ASSIGN dgemm:temp (* dgemm:alpha dgemm:a(dgemm:i , dgemm:l))
      IF (== dgemm:beta 1.0000000000000000e0_8)
        RETURN
      ENDIF
""")


# ---------------------------------------------------------------------------
# Synthetic dump for a COMPLEX*16 subroutine (like ZGEMM)
# ---------------------------------------------------------------------------

ZGEMM_DUMP = textwrap.dedent("""\
    Namespace: A-H: (REAL 4) I-N: (INTEGER 4) O-Z: (REAL 4)
    procedure name = zgemm
      symtree: 'a'           || symbol: 'a'
        type spec : (COMPLEX 8)
        attributes: (VARIABLE  DIMENSION DUMMY)
      symtree: 'alpha'       || symbol: 'alpha'
        type spec : (COMPLEX 8)
        attributes: (VARIABLE  DUMMY)
      symtree: 'dconjg'      || symbol: 'dconjg'
        type spec : (COMPLEX 8)
        attributes: (PROCEDURE INTRINSIC-PROC  INTRINSIC FUNCTION ELEMENTAL PURE)
        result: dconjg
      symtree: 'zgemm'       || symbol: 'zgemm'
        type spec : (UNKNOWN 0)
        attributes: (PROCEDURE  SUBROUTINE ARRAY-OUTER-DEPENDENCY)
        Formal arglist: transa transb m n k alpha a lda b ldb beta c ldc
      symtree: 'xerbla'      || symbol: 'xerbla'
        type spec : (UNKNOWN 0)
        attributes: (PROCEDURE EXTERNAL-PROC  EXTERNAL SUBROUTINE ARRAY-OUTER-DEPENDENCY)

      code:
      CALL xerbla (('ZGEMM ') (zgemm:info))
      ASSIGN zgemm:temp dconjg[[((zgemm:a(zgemm:l , zgemm:i)))]]
""")


# ---------------------------------------------------------------------------
# Synthetic dump for a FUNCTION (like DNRM2)
# ---------------------------------------------------------------------------

DNRM2_DUMP = textwrap.dedent("""\
    Namespace: A-H: (REAL 4) I-N: (INTEGER 4) O-Z: (REAL 4)
    procedure name = dnrm2
      symtree: 'abs'         || symbol: 'abs'
        type spec : (REAL 4)
        attributes: (PROCEDURE INTRINSIC-PROC  FUNCTION IMPLICIT-TYPE ARRAY-OUTER-DEPENDENCY)
        result: abs
      symtree: 'ax'          || symbol: 'ax'
        type spec : (REAL 8)
        attributes: (VARIABLE )
      symtree: 'dnrm2'       || symbol: 'dnrm2'
        type spec : (REAL 8)
        attributes: (PROCEDURE  FUNCTION)
        result: dnrm2
        Formal arglist: n x incx
      symtree: 'n'           || symbol: 'n'
        type spec : (INTEGER 4)
        attributes: (VARIABLE  DUMMY)
      symtree: 'sqrt'        || symbol: 'sqrt'
        type spec : (REAL 4)
        attributes: (PROCEDURE INTRINSIC-PROC  INTRINSIC FUNCTION IMPLICIT-TYPE ELEMENTAL PURE)
        result: sqrt

      code:
      ASSIGN dnrm2:ax (* 2.5000000000000000e-1_8 dnrm2:scl)
      ASSIGN dnrm2:dnrm2 __builtin_sqrt[[((dnrm2:ssq))]]
""")


# ===================================================================
# Tests
# ===================================================================


class TestParseTreeFactsSubroutine:
    """Tests using the DGEMM-like subroutine dump."""

    def setup_method(self):
        self.facts = parse_tree_facts(DGEMM_DUMP)

    def test_routine_def(self):
        assert len(self.facts.routine_defs) == 1
        rd = self.facts.routine_defs[0]
        assert rd.kind == 'subroutine'
        assert rd.name == 'DGEMM'
        assert rd.return_type is None

    def test_type_decls(self):
        dp_names = {td.names[0] for td in self.facts.type_decls
                    if td.type_spec == 'DoublePrecision'}
        assert 'A' in dp_names
        assert 'ALPHA' in dp_names
        assert 'TEMP' in dp_names

    def test_no_parameter_in_type_decls(self):
        """PARAMETER symbols (ONE, ZERO) should NOT appear as type decls."""
        all_names = {td.names[0] for td in self.facts.type_decls}
        assert 'ONE' not in all_names

    def test_call_sites(self):
        call_names = {cs.name for cs in self.facts.call_sites if cs.is_call_stmt}
        assert 'XERBLA' in call_names

    def test_function_references(self):
        ref_names = {cs.name for cs in self.facts.call_sites if not cs.is_call_stmt}
        assert 'LSAME' in ref_names

    def test_external_names(self):
        assert 'LSAME' in self.facts.external_names
        assert 'XERBLA' in self.facts.external_names

    def test_intrinsic_names(self):
        assert 'MAX' in self.facts.intrinsic_names

    def test_char_literals(self):
        assert 'DGEMM ' in self.facts.char_literals
        assert 'N' in self.facts.char_literals

    def test_real_literals(self):
        assert any('1.0000000000000000e0_8' in lit
                    for lit in self.facts.real_literals)


class TestParseTreeFactsComplex:
    """Tests using the ZGEMM-like complex dump."""

    def setup_method(self):
        self.facts = parse_tree_facts(ZGEMM_DUMP)

    def test_complex_type_decls(self):
        complex_decls = [td for td in self.facts.type_decls
                         if td.type_spec == 'Complex']
        assert len(complex_decls) >= 2
        assert all(td.kind_value == '16' for td in complex_decls)

    def test_intrinsic_dconjg(self):
        assert 'DCONJG' in self.facts.intrinsic_names

    def test_char_literal_zgemm(self):
        assert 'ZGEMM ' in self.facts.char_literals


class TestParseTreeFactsFunction:
    """Tests using the DNRM2-like function dump."""

    def setup_method(self):
        self.facts = parse_tree_facts(DNRM2_DUMP)

    def test_function_def(self):
        assert len(self.facts.routine_defs) == 1
        rd = self.facts.routine_defs[0]
        assert rd.kind == 'function'
        assert rd.name == 'DNRM2'
        assert rd.return_type == 'DoublePrecision'

    def test_variable_type_decl(self):
        names = {td.names[0] for td in self.facts.type_decls}
        assert 'AX' in names

    def test_intrinsics(self):
        assert 'SQRT' in self.facts.intrinsic_names

    def test_real_literals(self):
        assert any('2.5000000000000000e-1_8' in lit
                    for lit in self.facts.real_literals)


class TestInternalSymbolFiltering:
    """gfortran internal symbols (__max_i4, etc.) should be excluded."""

    def test_no_dunder_in_call_sites(self):
        facts = parse_tree_facts(DGEMM_DUMP)
        for cs in facts.call_sites:
            assert not cs.name.startswith('__'), \
                f"Internal symbol {cs.name} should be filtered"


class TestEmptyInput:
    def test_empty_string(self):
        facts = parse_tree_facts('')
        assert facts.routine_defs == []
        assert facts.type_decls == []
        assert facts.call_sites == []


# ---------------------------------------------------------------------------
# Integration test (requires gfortran on PATH)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(find_gfortran() is None,
                    reason='gfortran not found on PATH')
class TestGfortranIntegration:
    """Integration tests that invoke real gfortran."""

    BLAS_DIR = Path('external/lapack-3.12.1/BLAS/SRC')

    @pytest.mark.skipif(
        not Path('external/lapack-3.12.1/BLAS/SRC/dgemm.f').exists(),
        reason='BLAS source not available')
    def test_scan_dgemm(self):
        facts = scan_file(self.BLAS_DIR / 'dgemm.f')
        assert facts is not None
        assert any(rd.name == 'DGEMM' for rd in facts.routine_defs)
        assert any(td.type_spec == 'DoublePrecision'
                   for td in facts.type_decls)
        assert any(cs.name == 'XERBLA' for cs in facts.call_sites)
        assert 'LSAME' in facts.external_names

    @pytest.mark.skipif(
        not Path('external/lapack-3.12.1/BLAS/SRC/zgemm.f').exists(),
        reason='BLAS source not available')
    def test_scan_zgemm(self):
        facts = scan_file(self.BLAS_DIR / 'zgemm.f')
        assert facts is not None
        assert any(td.type_spec == 'Complex' for td in facts.type_decls)
        assert 'DCONJG' in facts.intrinsic_names

    @pytest.mark.skipif(
        not Path('external/lapack-3.12.1/BLAS/SRC/dnrm2.f90').exists(),
        reason='BLAS source not available')
    def test_scan_dnrm2(self):
        facts = scan_file(self.BLAS_DIR / 'dnrm2.f90')
        assert facts is not None
        assert any(rd.kind == 'function' and rd.name == 'DNRM2'
                   for rd in facts.routine_defs)
