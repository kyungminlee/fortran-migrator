"""Tests for ``route_sources`` — the common/precision archive router.

``cmd_stage`` (``staging._route_manifest_sources``) and ``cmd_build``
(``__main__._run_build``) once carried two near-identical copies of this
routing chain. They now share one implementation and differ only in four
keyword arguments, so what needs pinning is (a) the chain's priority
order and (b) each of those four deltas — a silent change to either
re-files sources between ``*_common`` and the precision archive, which
the byte-identity gate can only catch after the fact.

The chain, highest priority first: ``force_common`` → ``precision_stems``
→ ``copy_files`` (when routed to precision) → ``independent`` →
``copy_files`` (when routed to common) → precision. Foreign LA helper
pairs are dropped before any of it.
"""

from pathlib import Path
from types import SimpleNamespace

import pytest

from migrator.cli_common import route_sources
from migrator.staging import _route_manifest_sources
from migrator.target_mode import load_target

STAGE_REL = (lambda f: f'src/{f.name}')


def _config(force_common=(), copy_files=()):
    """Minimal stand-in for a RecipeConfig — the router reads two sets,
    both already upper-cased by ``config.load_recipe``."""
    return SimpleNamespace(force_common=set(force_common),
                           copy_files=set(copy_files))


def _paths(*names):
    return [Path('/out') / n for n in names]


def _route(files, *, target='kind16', force_common=(), copy_files=(),
           independent=(), **kwargs):
    kwargs.setdefault('rel', STAGE_REL)
    return route_sources(files, load_target(target),
                         _config(force_common, copy_files),
                         set(independent), **kwargs)


# --- the chain ------------------------------------------------------

def test_default_route_is_precision():
    common, precision = _route(_paths('dgemm.f'))
    assert common == []
    assert precision == ['src/dgemm.f']


def test_independent_routes_to_common():
    common, precision = _route(_paths('lsame.f'), independent={'LSAME'})
    assert common == ['src/lsame.f']
    assert precision == []


def test_force_common_beats_everything():
    """The recipe override outranks the scanner *and* the caller's own
    precision list — BLACS pins the integer entry points this way."""
    common, precision = _route(
        _paths('igamx2d.c'), force_common={'IGAMX2D'},
        copy_files={'IGAMX2D'}, precision_stems={'IGAMX2D'})
    assert common == ['src/igamx2d.c']
    assert precision == []


def test_precision_stems_beats_independent():
    """A rename target (``disnan.f`` → ``qisnan.f``) can collide with a
    stem the scanner called independent; it must still be precision."""
    common, precision = _route(_paths('qisnan.f'), independent={'QISNAN'},
                               precision_stems={'QISNAN'})
    assert common == []
    assert precision == ['src/qisnan.f']


def test_order_is_preserved_within_each_list():
    common, precision = _route(
        _paths('a.f', 'b.f', 'c.f', 'd.f'), independent={'A', 'C'})
    assert common == ['src/a.f', 'src/c.f']
    assert precision == ['src/b.f', 'src/d.f']


# --- delta 1: rel ---------------------------------------------------

def test_rel_controls_the_manifest_spelling():
    """``cmd_stage`` writes ``src/<name>``; ``cmd_build`` writes a path
    relative to its output dir. Neither manifest may pick up the other's
    spelling."""
    files = [Path('/out/sub/dgemm.f')]
    _, staged = route_sources(files, load_target('kind16'), _config(), set(),
                              rel=STAGE_REL)
    _, built = route_sources(files, load_target('kind16'), _config(), set(),
                             rel=lambda f: str(f.relative_to(Path('/out'))))
    assert staged == ['src/dgemm.f']
    assert built == ['sub/dgemm.f']


# --- delta 2: strip_trailing_underscore -----------------------------

def test_strip_trailing_underscore_matches_the_undecorated_name():
    """BLACS/SRC holds 78 ``*_.c`` stems; ``cmd_stage`` lets the recipe
    list the logical name, exactly as ``skip_files`` does."""
    common, precision = _route(_paths('igamx2d_.c'),
                               force_common={'IGAMX2D'},
                               strip_trailing_underscore=True)
    assert common == ['src/igamx2d_.c']
    assert precision == []


def test_without_the_strip_the_decorated_stem_does_not_match():
    """``cmd_build`` must not inherit it, or its ``force_common``
    entries would re-route."""
    common, precision = _route(_paths('igamx2d_.c'),
                               force_common={'IGAMX2D'})
    assert common == []
    assert precision == ['src/igamx2d_.c']


def test_strip_still_matches_the_decorated_spelling():
    common, _ = _route(_paths('igamx2d_.c'), force_common={'IGAMX2D_'},
                       strip_trailing_underscore=True)
    assert common == ['src/igamx2d_.c']


# --- delta 3: copy_files_to -----------------------------------------

def test_copy_files_to_precision_beats_independent():
    common, precision = _route(_paths('la_constants_qx.F90'),
                               copy_files={'LA_CONSTANTS_QX'},
                               independent={'LA_CONSTANTS_QX'},
                               copy_files_to='precision')
    assert common == []
    assert precision == ['src/la_constants_qx.F90']


def test_copy_files_to_common_loses_to_independent():
    """The branch position, not just the destination, differs between
    the two callers — ``'common'`` sits *after* the scanner's verdict."""
    common, precision = _route(_paths('la_constants_qx.F90'),
                               copy_files={'LA_CONSTANTS_QX'},
                               independent={'LA_CONSTANTS_QX'},
                               copy_files_to='common')
    assert common == ['src/la_constants_qx.F90']
    assert precision == []


def test_copy_files_to_common_catches_stems_the_scanner_never_saw():
    """A Fortran ``copy_files`` entry in a C recipe never appears in
    ``independent``; without the explicit route it would fall through to
    precision."""
    common, precision = _route(_paths('pilaenv.f'), copy_files={'PILAENV'},
                               copy_files_to='common')
    assert common == ['src/pilaenv.f']
    assert precision == []


def test_copy_files_to_rejects_anything_else():
    with pytest.raises(ValueError, match='copy_files_to'):
        _route(_paths('dgemm.f'), copy_files_to='both')


# --- foreign LA helper pairs ----------------------------------------

def test_foreign_la_pairs_are_dropped():
    """All three pairs share LAPACK's SRC dir; the ``*_mw`` pair does
    ``use multifloats``, which a kind16 build has no module for."""
    files = _paths('la_constants_qx.F90', 'la_constants_ey.F90',
                   'la_xisnan_mw.F90', 'dgemm.f')
    common, precision = _route(files, target='kind16',
                               independent={'LA_CONSTANTS_QX'})
    assert common == ['src/la_constants_qx.F90']
    assert precision == ['src/dgemm.f']


def test_which_pair_is_foreign_follows_the_target():
    files = _paths('la_constants_qx.F90', 'la_constants_ey.F90')
    _, kept_for_kind10 = _route(files, target='kind10')
    _, kept_for_kind16 = _route(files, target='kind16')
    assert kept_for_kind10 == ['src/la_constants_ey.F90']
    assert kept_for_kind16 == ['src/la_constants_qx.F90']


def test_force_common_cannot_resurrect_a_foreign_pair():
    common, precision = _route(_paths('la_constants_ey.F90'),
                               target='kind16',
                               force_common={'LA_CONSTANTS_EY'})
    assert common == []
    assert precision == []


# --- the staging caller ---------------------------------------------

def _classification(independent=(), rename_map=None):
    return SimpleNamespace(
        independent=set(independent),
        build_rename_map=lambda target_mode: dict(rename_map or {}))


def test_staging_caller_forces_its_own_la_pair_to_precision():
    """Single-precision by construction: this target's E/Y or Q/X block
    belongs in libqlapack, never in ``lapack_common`` — even though the
    scanner marks it independent and it is a ``copy_files`` entry."""
    common, precision = _route_manifest_sources(
        _paths('la_constants_qx.F90'),
        _config(copy_files={'LA_CONSTANTS_QX'}),
        load_target('kind16'),
        _classification(independent={'LA_CONSTANTS_QX'}))
    assert common == []
    assert precision == ['src/la_constants_qx.F90']


def test_staging_caller_forces_rename_targets_to_precision():
    common, precision = _route_manifest_sources(
        _paths('qisnan.f'),
        _config(),
        load_target('kind16'),
        _classification(independent={'QISNAN'},
                        rename_map={'DISNAN': 'qisnan'}))
    assert common == []
    assert precision == ['src/qisnan.f']


def test_staging_caller_sends_leftover_copy_files_to_common():
    common, precision = _route_manifest_sources(
        _paths('pilaenv.f'),
        _config(copy_files={'PILAENV'}),
        load_target('kind16'),
        _classification())
    assert common == ['src/pilaenv.f']
    assert precision == []


def test_staging_caller_strips_the_c_decoration_underscore():
    common, precision = _route_manifest_sources(
        _paths('igamx2d_.c'),
        _config(force_common={'IGAMX2D'}),
        load_target('kind16'),
        _classification())
    assert common == ['src/igamx2d_.c']
    assert precision == []
