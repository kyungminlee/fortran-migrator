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

## pdlanhs / pxlanhs — FIXED

- Root cause was an upstream ScaLAPACK 2.2.3 bug in pdlanhs.f /
  pzlanhs.f: the NPROW=1 path failed to advance ``II`` after the
  first column block, so the inner-loop bound ``MIN(II+LL-JJ+1,
  IIA+NP-1)`` collapsed to row 2 for every column past the first
  MB. M-norm passed by luck (max element typically lay in the kept
  range); 1/F/I norms underestimated by 10-20%. Fixed via
  ``source_overrides`` in recipes/scalapack.yaml — the patched
  upstream bodies live in recipes/scalapack/source_overrides/ and
  go through the normal migration pipeline.
- Test drivers: tests/scalapack/auxiliary/test_pdlanhs.f90 /
  test_pzlanhs.f90 (both PASS to full target precision on all 3
  targets).

## pxhetrd / pxheev — FIXED

Originally these returned garbage (rel err ~3-5) on the Hermitian
side while their real-symmetric counterparts (pdsytrd / pdsyev) passed
cleanly. Root cause was the precision-independent PBLAS dispatcher
files (PB_C* helpers) shipping with stale ``(double*)`` pointer-cast
strides on KIND targets — fixed by the std/extension archive split,
which routes the alias-widened dispatcher through the migrated
archive while the standard double-precision body lives separately.
Both ``test_pzhetrd`` (matrix-element residual) and ``test_pzheev``
(eigenvalue + residual on JOBZ='V') now pass to full target
precision (32+ digits on kind16).

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

(pzheev fix folded into the pxhetrd/pxheev entry above.)

## pdormrz / pzunmrz — descriptor-alignment failure

- **Symptom**: A test driver that invokes pdtzrzf followed by pdormrz('L', 'N', n_a, ncc, m_a, n_a-m_a, …) aborts in the internal PBETRAN dispatcher with "parameter number 11 had an illegal value", suggesting the m-by-n trapezoidal A's column distribution must align with C's row distribution in a way the simple "same MB/NB" descriptor pair doesn't satisfy.
- **Action**: Defer until the descA/descC alignment study is done. The wrappers are exposed (target_pdormrz / target_pzunmrz) and compile, so re-enabling is just a driver-shape fix.
