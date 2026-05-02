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

2×3 / 3×2 grids with `MIN(LCMQ,CEIL(M,NB)) > 1` covered 2026-05-02
via `test_pbdtran_lcm6.f90` (2×3, ADIST='C', M=12, N=NB=4 ⇒ LCMQ=2,
CEIL(M,NB)=3, MIN=2) and `test_pbztran_lcm6.f90` (3×2, ADIST='R',
M=NB=4, N=12 ⇒ LCMP=2, CEIL(N,NB)=3, MIN=2; TRANS='C'). These need
a connected 6-rank MPI world, which the default sandbox can't
provide; the harness exposes a separate `PBBLAS_TEST_NPROC_LCM6=6`
knob for filenames matching `*_lcm6` and the tests self-skip with a
"skipped" pseudo-case otherwise. Verified on the linux-impi preset
(Intel MPI with FI_PROVIDER=tcp + I_MPI_HYDRA_BOOTSTRAP=fork) for
all three targets — kind16 and kind10 (max_rel_err ≤ 3.3e-20) and
multifloats (max_rel_err ≤ 6.6e-33). These exercise the LCM-cycle
DO loops in `pbdtran.f` (DO 20, DO 30) and the `LCMQ.GT.1` /
`PBDTRGET` / `PBDTRSRT` paths that a 2×2 grid (LCMP=LCMQ=1)
collapses past.

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
