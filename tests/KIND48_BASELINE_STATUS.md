# kind4 / kind8 baseline columns — status

## Where we are

The structural pieces for producing `kind4` and `kind8` baseline columns
in `tests/RESULT.md` are wired up:

- `scripts/precision_report.py:30` — `TARGET_ORDER` lists kind4 / kind8
  first, so when JSON reports are stitched the columns appear before
  kind10.
- `tests/common/target_kind{4,8}/target_conv.f90` — quad↔native-precision
  converters (real(4)↔real(16) and real(8)↔real(16)).
- `tests/{blas,lapack,pblas,scalapack}/target_kind{4,8}/target_*.fypp` —
  per-framework wrapper shims setting prefix `s/c` (kind4) and `d/z`
  (kind8); pblas/scalapack use `ps/pc` and `pd/pz`.
- `targets/kind{4,8}.yaml` — TargetMode YAMLs (skeleton; not actually
  consumed because cmd_stage short-circuits).
- `src/migrator/__main__.py:_stage_baseline` — `migrator stage --target
  kind4|kind8` skips migration entirely and writes a target_config.cmake
  with `LIB_PREFIX=""`.
- `cmake/CMakeLists.txt` — `BASELINE_BUILD` flag derived from empty
  LIB_PREFIX; every `add_migrated_*` call and the install loop are gated
  on `NOT BASELINE_BUILD`.
- `tests/{blas,lapack,pblas,scalapack}/CMakeLists.txt` and
  `tests/CMakeLists.txt` — under BASELINE_BUILD the test wrappers link
  the bare standard archives (`blas`, `lapack`, `pblas`, `scalapack`)
  instead of `${LIB_PREFIX}<lib>`.

`migrator stage --target kind8 /tmp/x && cmake -S /tmp/x -B /tmp/x/build &&
cmake --build /tmp/x/build --target test_dasum -j8` succeeds. The
executable links and runs.

## Why it doesn't yet produce valid numbers

`librefblas_quad.a` (compiled with `-freal-8-real-16`) and the std
`libblas-gfortran-13.a` BOTH define `dasum_`, `dgemm_`, etc. — same
symbol, different ABI (one takes/returns REAL(KIND=16), the other
REAL(KIND=8)).

For migrated targets (kind10/16/multifloats) this is hidden by link
order: `refblas_quad` is named before the std archive on the test
executable's link line, so `dasum_` resolves from refblas_quad and the
std archive's `dasum_` is never extracted. The migrated wrapper calls
`qasum_` / `easum_` / etc., which don't collide with anything.

For the kind4/kind8 baseline the wrapper calls the bare native symbol
(`dasum_` for kind8, `sasum_` for kind4) — same name as the quad
reference. Whichever archive's `dasum_` ld picks gets called by *both*
the quad-reference path (via `tests/blas/common/ref_quad_blas.f90`'s
interface declarations, which expect REAL(KIND=16) ABI) and the
under-test path (via `target_blas_body.fypp`'s interface declarations,
which expect REAL(KIND=8) ABI). One of them is ABI-mismatched →
SIGSEGV and garbage output.

Confirmed: smoke run of `blas_test_dasum` in a kind8 baseline build
produces `max_rel_err ≈ 1.4e+169` and segfaults.

## What's needed to finish

The reference archives need to expose their quad-promoted symbols under
distinct names so the std archive's native-precision symbols can coexist
in the same link. Concretely:

1. `tests/blas/refblas/CMakeLists.txt` and
   `tests/lapack/reflapack/CMakeLists.txt` — post-build step (e.g.
   `objcopy --redefine-syms`) renames every Fortran-mangled global
   symbol in `librefblas_quad.a` / `libreflapack_quad.a` from `<name>_`
   to `<name>_quad_`.
2. `tests/blas/common/ref_quad_blas.f90`,
   `tests/lapack/common/ref_quad_lapack.f90`,
   `tests/pblas/common/ref_quad_blas.f90`,
   `tests/scalapack/common/ref_quad_{blas,lapack}.f90` — every
   `function`/`subroutine` interface declaration gains
   `bind(c, name='<name>_quad')` so calls go to the renamed symbols.
   ~750 routine declarations across ~9400 lines.
3. Verify migrated tests (kind10/kind16/multifloats) still pass — they
   don't care about the rename mechanically, but the mass-edit could
   trip on edge cases (character arguments need
   `character(kind=c_char)`, etc.).
4. Stage and build kind4 + kind8 trees, run tests, regenerate
   `tests/RESULT.md` with five columns.

Step 2 is the bulk of the work and benefits from a careful scripted
edit rather than hand-by-hand. Suggested approach: a small Python
script that walks each ref_quad_*.f90, finds every `function NAME(...)`
or `subroutine NAME(...)` line at column 9 (indented inside an
interface block), inserts `bind(c, name='NAME_quad')` before any
trailing `result(...)`. Run the migrated tests as a regression check
after the edit.

