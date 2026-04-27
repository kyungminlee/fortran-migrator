# tests/pbblas — Outstanding Work

## Untested PBBLAS routines

These ship in `external/scalapack-2.2.3/PBLAS/SRC/PBBLAS/` and are
exercised at runtime only as part of the wider PBLAS pipeline; we do
not cover them here.

### `pbdtrnv` / `pbztrnv` — distributed vector transpose

Exposes the same column↔row reallocation as `pbdtran` but for
vectors. The interface takes an `NZ` block-offset argument, two
broadcast-control row/col pairs (`IXROW/IXCOL`, `IYROW/IYCOL`), and a
`WORK` buffer whose minimum size is documented only obliquely in the
upstream comment block. A wrong test for this routine is worse than no
test, so it's deferred until either:

  - we copy the PBLAS-level call site that drives it (e.g. `PB_Cpgemv`)
    and use it as a reference fixture, or
  - someone writes the Mqb/Npb/LCM/IGD bookkeeping by hand and verifies
    it against ScaLAPACK's own test driver.

The wrapper interfaces (`target_pbdtrnv` / `target_pbztrnv`) are
already wired in `common/target_pbblas_body.fypp` so a future test
program can `use target_pbblas` and call them directly.

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
  call the `tr*` helpers. We supplied byte-for-byte re-implementations
  in `common/scalapack_tools.f90` so the pbblas test executables
  resolve at link time.

  Cleaner fix (out of scope for this subtree):
  - either add a `scalapack_tools` recipe (just the three integer
    helpers) and depend on it from pbblas / pblas, or
  - declare the dependency in `recipes/pbblas.yaml` so the migrator
    pulls them into libqpbblas / libqpblas at staging time.

  Until then, the duplicated source in `common/scalapack_tools.f90`
  remains the simplest path that respects the "tests-only" constraint.

## Multi-target coverage

This subtree only ships a `target_kind16/` shim. `target_kind10/` and
`target_multifloats/` should be straightforward copies once the kind10
and multifloats stagings are exercised — model on tests/ptzblas which
already supports all three.
