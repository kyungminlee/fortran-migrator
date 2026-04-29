# tests/scalapack — known upstream / migrator gaps

## Banded solver families (packed-banded layout helper)

- **Status**: Tridiagonal piece landed. The 1D BLACS context and
  `descinit_1d` helper were added to `pblas_grid.f90`; tests cover
  `pddtsv` and `pdptsv` (general + PD tridiagonal). The remaining
  banded entries (`pdgbsv`/`pdgbtrf`/`pdgbtrs`, `pdpbsv`/`pdpbtrf`/
  `pdpbtrs`, `pddbsv`/`pddbtrf`/`pddbtrs`, plus complex mirrors) are
  still deferred.
- **What's needed for banded**: ScaLAPACK's banded format stores `A`
  in `(2*bwl+2*bwu+1) × LOCc` packed strips with extra rows reserved
  for fill-in during D&C factorisation, which differs from the LAPACK
  banded layout. A small layout helper that packs from the natural
  diagonal-by-diagonal representation into ScaLAPACK's storage (and
  unpacks the reference LAPACK layout for `dgbsv`/`dpbsv` comparison)
  would unblock those tests.
- **Action**: Implement the packed-banded helper as a sibling to
  `descinit_1d`, then add per-family drivers (general → PD → DD
  banded). Each pulls in the same 1D scaffold the tridiagonal tests
  use.



Routines whose differential-precision test would otherwise live in this
suite but cannot pass with the current toolchain. Each entry documents
the symptom so a future pass can re-enable a test once the underlying
issue is fixed outside `tests/scalapack/`.

## pdsyevx / pzheevx / pdsygvx / pzhegvx (RANGE-selector eigensolvers)

- **Status**: The pdlaiect symbol clash that originally blocked these
  tests has been dissolved by the std/extension archive split (see
  the commit history around recipes/scalapack_c.yaml's `extra_renames`
  field). Standard archive `scalapack_c` now owns the upstream
  `pdlaiectb_/pdlaiectl_`; migrated archive `${LIB_PREFIX}scalapack_c`
  owns the renamed `pqlaiectb_/pqlaiectl_` (etc.). The wrappers
  `target_pdsyevx`/`target_pzheevx`/`target_pdstebz` are now in the
  shared template (`tests/scalapack/common/target_scalapack_body.fypp`).
- **Remaining symptom**: A `test_pdsyevx` driver computes eigenvalues
  to full target precision (n=32 reports max_rel_err ≈ 1e-33 on
  kind16, ~1e-19 on kind10, ~4e-32 on multifloats — all PASS the
  precision check) but every PDSYEVX call leaves the heap in a
  corrupted state — `free(): invalid pointer` aborts on first
  deallocate after the call. The corruption is per-call (does not
  require multiple sizes to surface) and target-independent. PDSYEV /
  PDSYEVR / PDSYEVD on the same scaffold pass cleanly on all sizes,
  so the issue is specific to PDSYEVX's bisection + inverse-iteration
  workspace handling.
- **Action**: Either find the upstream bug (likely an out-of-bounds
  write in PDSYEVX's eigenvalue-cluster bookkeeping that corrupts
  glibc's malloc free-lists) or write a custom safety harness that
  catches/sandboxes the corruption so the test can still report
  precision results. PDSTEBZ alone (the eigenvalue-only piece, no
  eigenvectors) might be testable in isolation without the broken
  cleanup.

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

## pdgecon / pdpocon (condition-number estimators)

- **Symptom**: `target_pdgecon` / `target_pdpocon` wrappers run cleanly
  but a straight `rel_err_scalar(rcond_scalapack, rcond_lapack)`
  comparison drifts to 50-200% disagreement, scaling roughly with n.
- **Diagnosis**: ScaLAPACK's `PDLACON` and LAPACK's `DLACON` are
  iterative 1-norm estimators (Hager-Higham). Same input matrix,
  different parallel ordering / starting vectors → different estimates.
  This is algorithmic non-determinism, not a precision bug, so a
  `target_eps`-driven tolerance does not apply.
- **Action**: A meaningful test would compute true `kappa_1(A)` from a
  quad SVD of A on rank 0, then verify the documented LACON guarantee
  `kappa_true / 3 <= kappa_est <= kappa_true`. Deferred until that
  helper exists; the wrappers remain exposed.

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
