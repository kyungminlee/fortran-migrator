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
tests run on. Coverage gaps still open:

  - NZ > 0 (block-offset start) — the upstream NZ bookkeeping in
    pbdtrnv.f is not exercised.
  - XDIST='R' (row-vector input, column-vector output).
  - IXCOL=-1 / IYROW=-1 replicated paths.
  - Larger / non-square grids where `LCMP > 1` or `LCMQ > 1` (the
    sandbox's mpiexec produces unconnected MPI worlds, so the
    distributed paths degenerate to local).

Closing these would benefit from copying the PBLAS-level call site
that drives pbdtrnv (e.g. `PB_Cpgemv`) and using it as a reference
fixture rather than open-coding the LCM/IGD bookkeeping.

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

## Multi-target coverage

This subtree only ships a `target_kind16/` shim. `target_kind10/` and
`target_multifloats/` should be straightforward copies once the kind10
and multifloats stagings are exercised — model on tests/ptzblas which
already supports all three.
