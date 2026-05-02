# tests/pbblas — Outstanding Work

## Untested PBBLAS routines

These ship in `external/scalapack-2.2.3/PBLAS/SRC/PBBLAS/` and are
exercised at runtime only as part of the wider PBLAS pipeline; we do
not cover them here.

### `pbdtrnv` / `pbztrnv` — covered

`test_pbdtrnv.f90` / `test_pbztrnv.f90` (baseline XDIST='C', NZ=0,
2×2 grid) plus four families of follow-ups, all added 2026-05-02
and passing on every target (kind10 / kind16 / multifloats):

  - `test_pb[dz]trnv_xdistr.f90` — XDIST='R' (row-vector input →
    column-vector output) on the default 2×2 grid.
  - `test_pb[dz]trnv_replicated.f90` — IXCOL=-1 and IYROW=-1
    simultaneously; X populated on every process column, Y verified
    on row 0 against the reference and cross-checked byte-equal
    against row `nprow-1`.
  - `test_pb[dz]trnv_nzoffset.f90` — NZ > 0 (NB=4, NZ=2, N=12).
    Source-row / source-column local layouts honor the `numroc(NN)
    - NZ` shrink for the offset extended NN=N+NZ space; non-source
    processes hold full NB-sized blocks.
  - `test_pb[dz]trnv_lcm.f90` — non-square grid via the new
    `grid_init_shape` helper. `pbdtrnv_lcm` runs 4×1 (LCMQ=4) with
    XDIST='C'; `pbztrnv_lcm` runs 1×4 (LCMP=4) with XDIST='R'. Both
    require a connected 4-rank MPI world: skipped (with a "skipped"
    pseudo-case) when the launcher spawns unconnected worlds.
    Verified on the linux-impi preset with `FI_PROVIDER=tcp` +
    `I_MPI_HYDRA_BOOTSTRAP=fork`.

No remaining coverage gaps for `pbdtrnv` / `pbztrnv`.

### `pbdtran` / `pbztran` — covered

Baseline `test_pbdtran.f90` / `test_pbztran.f90` cover ADIST='C'
with IACOL=0 and ADIST='R' with IAROW=0. Replicated-A added
2026-05-02 via `test_pb[dz]tran_replicated.f90` — both ADIST='C'
with IACOL=-1 (every process column holds A) and ADIST='R' with
IAROW=-1 (every process row holds A). These exercise the separate
"all column processors have a copy" branch in pbdtran.f (line ~352
onward) which uses `PBDTRGET` / `PBDTRSRT` rather than the
source-node-sends pattern. All three targets pass on both the
sandbox 2×2 launch and a real connected 2×2 via the linux-impi
preset.

The inner LCM-cycle logic still collapses on a 2×2 grid (LCMP=
LCMQ=1); the WORK-formula factor `MIN(LCMQ,CEIL(M,NB))` only goes
above 1 on 2×3 / 3×2 / 6+-rank grids, which the local
`PBBLAS_TEST_NPROC=4` harness can't produce. The branches
themselves are exercised — that was the gap the original TODO
bullet referred to.

### `pbdtrget` / `pbdtrsrt` / `pbdtrst1` (and complex variants)

Internal helpers used by `pbdtran` / `pbdtrnv` to gather, sort, and
stitch the diagonal-block representation of a column-block matrix.
Their inputs are intermediate buffers laid out by the caller and have
no standalone semantic meaning — testing them in isolation would
require fabricating a contrived data layout that may not match what
the caller actually produces. These are exercised transitively
through the now-comprehensive `pb[dz]tran` / `pb[dz]trnv` test
families (baseline + replicated + non-square-grid + NZ-offset), and
are best left at that coverage level until a regression motivates
digging in.

## Cross-tree follow-ups

- `external/scalapack-2.2.3/TOOLS/{numroc,iceil,ilcm}.f` are referenced
  by `pbdtran` / `pbztran` / `pbdtrnv` / `pbztrnv` but are not part of
  any of the five staged libraries (blas / blacs / ptzblas / pbblas /
  pblas). The PBLAS test programs link against pbblas only because the
  pbblas routines they happen to exercise (axpy, dot, ...) don't
  call the `tr*` helpers. We supplied behaviorally-equivalent
  re-implementations in `common/scalapack_tools.f90` so the pbblas
  test executables resolve at link time (see that file's header for
  the documented ICEIL/ILCM divergences from upstream).

  Cleaner fix (out of scope for this subtree):
  - either add a `scalapack_tools` recipe (just the three integer
    helpers) and depend on it from pbblas / pblas, or
  - declare the dependency in `recipes/pbblas.yaml` so the migrator
    pulls them into libqpbblas / libqpblas at staging time.

  Until then, the duplicated source in `common/scalapack_tools.f90`
  remains the simplest path that respects the "tests-only" constraint.
  Note: those reimplementations are *behaviorally equivalent on the
  documented input ranges* PBBLAS uses but are not line-for-line
  copies. NUMROC matches upstream; ICEIL adds an `IDENOM == 0 -> 0`
  guard; ILCM uses a `DO WHILE` Euclidean GCD rewrite and adds an
  `A == 0 .OR. B == 0 -> 0` guard. PBBLAS never calls these helpers
  with the guarded inputs, so observable behavior is identical for
  all callers in this subtree; the divergence is documented in the
  shim's header.

## Multi-target coverage — RESOLVED

All three target shims (`target_kind10/`, `target_kind16/`,
`target_multifloats/`) ship in this subtree and build cleanly via
`pbblas_test_target` on the corresponding stagings. Verified
2026-04-30. Runtime exercise still requires real MPI (the sandbox's
mpiexec spawns N unconnected worlds — see `tests/blacs/TODO.md`).
