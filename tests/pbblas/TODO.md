# tests/pbblas — Outstanding Work

## Untested PBBLAS routines

These ship in `external/scalapack-2.2.3/PBLAS/SRC/PBBLAS/` and are
exercised at runtime only as part of the wider PBLAS pipeline; we do
not cover them here.

### `pbdtrnv` / `pbztrnv` — covered, smoke-test only

`test_pbdtrnv.f90` and `test_pbztrnv.f90` exercise XDIST='C',
TRANS='T' (and TRANS='C' for the complex variant) with NZ=0 and a
generously-sized WORK buffer (`4*(N+NB)`) that envelopes the
upstream `CEIL(Nqb,LCMQ)*NB` formula on the 1×1 / 2×2 grids these
tests run on. XDIST='R' added 2026-05-02 via `test_pb[dz]trnv_xdistr.f90`
(row-vector input → column-vector output, 2×2 grid, all three
targets pass). Replicated paths (IXCOL=-1 and IYROW=-1 simultaneously)
added 2026-05-02 via `test_pb[dz]trnv_replicated.f90` — X populated
on every process column, Y verified on row 0 against the reference
and cross-checked byte-equal against row `nprow-1`; all three targets
pass. NZ > 0 added 2026-05-02 via `test_pb[dz]trnv_nzoffset.f90`
(NB=4, NZ=2, N=12; source-row / source-column local layouts honor
the `numroc(NN)-NZ` shrink for the offset extended NN=N+NZ space);
all three targets pass. Coverage gaps still open:

  - Larger / non-square grids where `LCMP > 1` or `LCMQ > 1` (the
    sandbox's mpiexec produces unconnected MPI worlds, so the
    distributed paths degenerate to local).

### `pbdtran` / `pbztran` — replicated-A (IACOL=-1) variant

`test_pbdtran.f90` covers ADIST='C' with IACOL=0 and ADIST='R' with
IAROW=0. The replicated paths (IACOL=-1 for ADIST='C', or IAROW=-1
for ADIST='R') require a different WORK size:
  Size(WORK) = N * CEIL(Mqb,LCMQ)*NB * MIN(LCMQ,CEIL(M,NB))
and a different in-WORK assembly via `PBDTRGET` / `PBDTRSRT`. These
paths only matter when every process column already holds a copy of
A; in 1×1 / 2×2 grids LCMQ=1 and the path collapses, so the added
test would not exercise its distinguishing logic on the sandbox.
Skipping until either a real distributed harness is available or a
regression motivates digging in.

### `pbdtrget` / `pbdtrsrt` / `pbdtrst1` (and complex variants)

Internal helpers used by `pbdtran` / `pbdtrnv` to gather, sort, and
stitch the diagonal-block representation of a column-block matrix.
Their inputs are intermediate buffers laid out by the caller and have
no standalone semantic meaning — testing them in isolation would
require fabricating a contrived data layout that may not match what
the caller actually produces. These are exercised transitively
through `test_pbdtran` / `test_pbztran` and are best left at that
coverage level until a regression motivates digging in.

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
