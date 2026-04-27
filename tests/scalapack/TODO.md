# tests/scalapack — known upstream / migrator gaps

Routines whose differential-precision test would otherwise live in this
suite but cannot pass with the current toolchain. Each entry documents
the symptom so a future pass can re-enable a test once the underlying
issue is fixed outside `tests/scalapack/`.

## pdsyevx / pzheevx (expert eigensolvers, RANGE selectors)

- **Symptom**: Linking a `test_pdsyevx` driver fails with multiple
  definitions of `pdlaiectb_` / `pdlaiectl_` between the migrated
  `libqscalapack_c.a` archive members `pqlaiect.c.o` and
  `pdlaiect.c.o`. Both compile units export the same C symbols even
  though one is the migrated rename of the other.
- **Diagnosis**: `pdlaiect.c` is a small C helper used by the
  expert symmetric eigensolver paths. The migrator's prefix
  classifier renames the file (`pq*`) but does not rename the
  exported function symbols, so the un-migrated `pd*` and the
  migrated `pq*` versions collide in the same archive.
- **Action**: Re-enable `tests/scalapack/eigenvalue/test_pdsyevx.f90`
  / `test_pzheevx.f90` once the migrator either drops the
  un-migrated `pdlaiect.c.o` from `libqscalapack_c` or renames its
  exported symbols. The wrappers (`target_pdsyevx`, `target_pzheevx`)
  also have to come back: they were initially landed alongside the
  test drivers, but the symbol clash poisons the link of *every*
  test that pulls in `scalapack_test_target`, so they had to be
  reverted out of the shared template.

## pdlanhs / pxlanhs (Hessenberg matrix norms)

- **Symptom**: For NORM='1' / 'F' / 'M' on n in {32,64,96}, the
  migrated `pdlanhs` / `pxlanhs` disagree with quad `dlanhs` /
  `zlanhs` at order 0.1 (only NORM='I' matches). Even after
  zeroing the input matrix below the first subdiagonal so the
  matrix is genuinely Hessenberg, the disagreement persists.
- **Diagnosis**: Likely a scaling or norm-of-norms bug specific to
  the migrated `pdlanhs` family — `pdlange`, `pdlansy`, `pdlantr`,
  `pzlange`, `pzlanhe`, `pzlantr` all agree to >30 digits with the
  same scaffold, so the harness is fine.
- **Action**: Re-enable `tests/scalapack/auxiliary/test_pdlanhs.f90`
  / `test_pzlanhs.f90` once the migrated `lanhs` family is fixed.
  Wrappers `target_pdlanhs` / `target_pzlanhs` remain in the
  template for later re-instatement.

## pxhetrd (complex Hermitian tridiagonal reduction)

- **Symptom**: A `test_pzhetrd` driver returns matrix-element errors of
  order 1 (`max_rel_err ~ 2.7-4.5`) on n in {32,64,96}, 2x2 grid,
  mb=nb=8 — totally wrong output. The real-symmetric counterpart
  `pdsytrd` passes to >32 digits on the same scaffold.
- **Diagnosis**: `pxhetrd` is the Hermitian analog of `pdsytrd` and
  is the back-end for `pxheev`. The same defect that makes `pxheev`
  return garbage eigenvalues likely lives here at the reduction
  level.
- **Action**: Re-enable `tests/scalapack/factorization/test_pzhetrd.f90`
  (wrapper `target_pzhetrd` is already in the template) when the
  underlying `pxhetrd` Householder path is fixed upstream.

## pxheev (complex Hermitian eigensolver, JOBZ='N' or 'V')

- **Symptom**: A test driver of the same shape as `test_pdsyev.f90`,
  hermitizing a random complex matrix and comparing eigenvalues
  against quad `zheev`, sees `max_rel_err ≈ 1–5` on eigenvalues and
  residual `||A·Z − Z·diag(W)||/||A|| ≈ 7–16` for both `JOBZ='N'`
  and `JOBZ='V'` across n ∈ {32, 64, 96}, 2×2 grid, mb=nb=8.
- **Diagnosis**: `target_pzheev` (the wrapper) is faithful to the
  ScaLAPACK 2.2 `PZHEEV` interface — `LWORK`/`LRWORK` are queried
  and obeyed; the input is correctly Hermitian. `pdsyev` with the
  identical scaffold passes to >32 digits, ruling out the harness.
  The eigenvalues returned by migrated `pxheev` simply do not agree
  with reference `zheev` on the same matrix.
- **Action**: Re-enable `tests/scalapack/eigenvalue/test_pzheev.f90`
  once the migrated `pxheev` is verified or replaced with the
  divide-and-conquer driver (`pxheevd`, also exposed but untested).
  Do not edit migrator/ from within this suite.
