# tests/scalapack — known upstream / migrator gaps

Routines whose differential-precision test would otherwise live in this
suite but cannot pass with the current toolchain. Each entry documents
the symptom so a future pass can re-enable a test once the underlying
issue is fixed outside `tests/scalapack/`.

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
