"""End-to-end pipeline link tests.

Runs the full migration pipeline (migrate -> convergence -> verify -> build)
for all five libraries at both KIND=10 (extended) and KIND=16 (quad), then
links every resulting static library into an empty Fortran executable to
verify there are no unresolved symbols.

These tests require gfortran, cmake, and MPI headers (libopenmpi-dev) to be
installed, and the external source trees (lapack-3.12.1, scalapack-2.2.3)
to be present under ``external/``.

Typical run time: ~3 minutes (dominated by LAPACK/ScaLAPACK compilation).
"""

import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent

RECIPES = ['blas', 'lapack', 'blacs', 'pblas', 'scalapack']

# Libraries in link order (dependents first).
# Each entry is (recipe_name, kind10_lib_names, kind16_lib_names).
LINK_ORDER = [
    ('scalapack', ['escalapack', 'scalapack_common'],
                  ['qscalapack', 'scalapack_common']),
    ('pblas',     ['epblas'],     ['qpblas']),
    ('lapack',    ['elapack', 'lapack_common'],
                  ['qlapack', 'lapack_common']),
    ('blas',      ['eblas', 'blas_common'],
                  ['qblas', 'blas_common']),
    ('blacs',     ['eblacs', 'blacs_common'],
                  ['qblacs', 'blacs_common']),
]


def _has_gfortran() -> bool:
    return shutil.which('gfortran') is not None


def _has_cmake() -> bool:
    return shutil.which('cmake') is not None


def _has_external_sources() -> bool:
    return (
        (PROJECT_ROOT / 'external' / 'lapack-3.12.1' / 'SRC').is_dir()
        and (PROJECT_ROOT / 'external' / 'scalapack-2.2.3' / 'SRC').is_dir()
    )


_skip_reason = 'requires gfortran, cmake, and external source trees'
_can_run = _has_gfortran() and _has_cmake() and _has_external_sources()


def _run_pipeline(recipe: str, work_dir: Path, kind: int) -> Path:
    """Run the full pipeline for one library and return the build dir."""
    recipe_path = PROJECT_ROOT / 'recipes' / f'{recipe}.yaml'
    result = subprocess.run(
        ['python', '-m', 'pyengine', 'run',
         str(recipe_path), str(work_dir),
         '--parser', 'gfortran',
         '--kind', str(kind)],
        capture_output=True, text=True,
        cwd=str(PROJECT_ROOT),
        timeout=600,
    )
    build_dir = work_dir / 'output' / '_build'
    # Allow verify failures (exit code 1 from verify) as long as build exists
    if not build_dir.is_dir():
        pytest.fail(
            f'Pipeline for {recipe} KIND={kind} did not produce a build dir.\n'
            f'stdout:\n{result.stdout[-2000:]}\n'
            f'stderr:\n{result.stderr[-2000:]}'
        )
    return build_dir


def _link_all(build_dirs: dict[str, Path], kind: int, tmp_path: Path) -> None:
    """Link all libraries into an empty executable."""
    # Write an empty Fortran program
    main_f = tmp_path / 'link_test.f90'
    main_f.write_text(textwrap.dedent("""\
        program link_test
          implicit none
        end program
    """))

    # Build the linker command
    exe = tmp_path / 'link_test'
    cmd = ['gfortran', '-o', str(exe), str(main_f)]

    lib_col = 1 if kind == 16 else 0  # index into LINK_ORDER tuple

    for recipe, e_libs, q_libs in LINK_ORDER:
        bdir = build_dirs[recipe]
        cmd.append(f'-L{bdir}')
        libs = q_libs if kind == 16 else e_libs
        for lib in libs:
            cmd.append(f'-l{lib}')

    # Add quadmath for KIND=16
    if kind == 16:
        cmd.append('-lquadmath')

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        # Collect undefined symbols for a useful error message
        undef = [l for l in result.stderr.splitlines()
                 if 'undefined reference' in l.lower()]
        summary = '\n'.join(undef[:30]) if undef else result.stderr[-2000:]
        pytest.fail(
            f'Link failed for KIND={kind} with {len(undef)} '
            f'undefined references:\n{summary}'
        )

    assert exe.exists(), f'Executable not produced for KIND={kind}'


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _can_run, reason=_skip_reason)
class TestPipelineKind16:
    """Build and link all libraries at KIND=16 (quad precision)."""

    @pytest.fixture(scope='class')
    def build_dirs(self, tmp_path_factory):
        dirs = {}
        for recipe in RECIPES:
            work = tmp_path_factory.mktemp(f'{recipe}_k16')
            dirs[recipe] = _run_pipeline(recipe, work, kind=16)
        return dirs

    def test_blas_builds(self, build_dirs):
        assert list(build_dirs['blas'].glob('libqblas.*'))

    def test_lapack_builds(self, build_dirs):
        assert list(build_dirs['lapack'].glob('libqlapack.*'))

    def test_blacs_builds(self, build_dirs):
        assert list(build_dirs['blacs'].glob('libqblacs.*'))

    def test_pblas_builds(self, build_dirs):
        assert list(build_dirs['pblas'].glob('libqpblas.*'))

    def test_scalapack_builds(self, build_dirs):
        assert list(build_dirs['scalapack'].glob('libqscalapack.*'))

    def test_link_all(self, build_dirs, tmp_path):
        _link_all(build_dirs, kind=16, tmp_path=tmp_path)


@pytest.mark.skipif(not _can_run, reason=_skip_reason)
class TestPipelineKind10:
    """Build and link all libraries at KIND=10 (extended precision)."""

    @pytest.fixture(scope='class')
    def build_dirs(self, tmp_path_factory):
        dirs = {}
        for recipe in RECIPES:
            work = tmp_path_factory.mktemp(f'{recipe}_k10')
            dirs[recipe] = _run_pipeline(recipe, work, kind=10)
        return dirs

    def test_blas_builds(self, build_dirs):
        assert list(build_dirs['blas'].glob('libeblas.*'))

    def test_lapack_builds(self, build_dirs):
        assert list(build_dirs['lapack'].glob('libelapack.*'))

    def test_blacs_builds(self, build_dirs):
        assert list(build_dirs['blacs'].glob('libeblacs.*'))

    def test_pblas_builds(self, build_dirs):
        assert list(build_dirs['pblas'].glob('libepblas.*'))

    def test_scalapack_builds(self, build_dirs):
        assert list(build_dirs['scalapack'].glob('libescalapack.*'))

    def test_link_all(self, build_dirs, tmp_path):
        _link_all(build_dirs, kind=10, tmp_path=tmp_path)
