# tests/scalapack — known upstream / migrator gaps

## Real-MPI exercise via `linux-impi` preset — 2026-04-30

`cmake --preset=linux-impi` invokes Intel MPI 2021.18 with
`I_MPI_ADJUST_REDUCE=1` injected via `MPIEXEC_PREFLAGS`. ScaLAPACK
tests went from `0% real signal` (singleton MPI) to `90%+ pass` on a
real 2×2 grid. 26 ScaLAPACK tests still fail (out of ~170):

```
banded   : pddbtrf pddbtrs pddbtrsv pdgbsv pdgbtrf pdgbtrs pdpttrf
           pzdbtrf pzdbtrs pzdbtrsv pzgbsv pzgbtrf pzgbtrs
refine   : pdgerfs pdporfs pdtrrfs   pzgerfs pzporfs pztrrfs
expert   : pdgesvx                   pzgesvx
norms    : pdlanhs                   pzlanhs
QR-fac   : pdggrqf                   pzggrqf
equiv    :                           pzgeequ
```

These all pass cleanly under impi's collectives once the reduce
shortpath is bypassed; the remaining failures look like banded-matrix
descriptor / block-size assumptions that don't hold under the real
distributed path, plus iterative-refinement convergence checks that
need a tolerance review now that the inner residual is computed via
real MPI rather than a singleton fallback.

Tracked here — these need per-routine investigation, not blanket
preset tweaks. For most of these, the fix is in the migrated
ScaLAPACK or the test reference, not in the harness.

## Banded solver families — DELIVERED

Tridiagonal piece (`pddtsv`/`pdptsv`) plus all 22 banded drivers
(`p[dz]gb{sv,trf,trs}`, `p[dz]pb{sv,trf,trs,trsv}`,
`p[dz]db{sv,trf,trs,trsv}`) landed in phases 45/46/55/56. Each driver
inlines its own packing — no shared `banded_pack` module materialized;
the per-driver loops are 3 lines and vary in `bw` vs `bwl/bwu` vs
`uplo`, so consolidation has no payoff. Diagonal row in PDGBTRF
storage: `BWL+2*BWU+1` (verified in `test_pdgbsv.f90`); PDDBTRF /
PDPBTRF use the standard LAPACK packing (diagonal at `BWU+1` /
`BW+1`).

Routines whose differential-precision test would otherwise live in this
suite but cannot pass with the current toolchain. Each entry documents
the symptom so a future pass can re-enable a test once the underlying
issue is fixed outside `tests/scalapack/`.

## pdsyevx / pzheevx — FIXED

Root cause was upstream PJLAENV's strict S/D/C/Z precision-letter gate
(external/scalapack-2.2.3/SRC/pjlaenv.f line 201-204) that returns the
function result *uninitialised* when the routine name doesn't start
with one of those letters. After migration, callers pass names like
'PESYTTRD' / 'PXHETTRD' / 'PMSYTTRD' — the gate trips, ANB picks up a
stack-resident garbage integer, NSYTRD_LWOPT overflows to INT32_MIN,
the caller allocates a tiny WORK and the next call scribbles past it.
Fixed via PJLAENV_EP (recipes/scalapack/extras/pjlaenv_ep.f) +
extra_renames PJLAENV → PJLAENV_EP in scalapack.yaml. The same
upstream bug affects LAPACK's ILAENV / ILAENV2STAGE / IPARAM2STAGE;
those got the parallel ILAENV_EP / ILAENV2STAGE_EP / IPARAM2STAGE_EP
treatment in lapack.yaml.

A second (independent) bug surfaced after PJLAENV_EP unblocked the
workspace query: target_pdsyevx and target_pzheevx allocated WORK_T(1)
(or RWORK_T(1)) for the LQUERY call, but upstream P{D,Z}SYEVX writes
WORK(1:3) on rank 0 (ABSTOL/VL/VU broadcast) before reporting the
optimal size. AddressSanitizer caught it; bumped to (R)WORK_T(3) in
the wrapper.

Test drivers: tests/scalapack/eigenvalue/test_pdsyevx.f90 and
test_pzheevx.f90 (RANGE='A', JOBZ='N'). Both PASS on kind10 / kind16 /
multifloats to ~18 digits.

pdsygvx / pzhegvx (the generalized eigensolver siblings) sit on the
same code path and are unblocked by the same fix; drivers TBD.

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

## *trsv (single-half banded/tridiag triangular solves) — driver semantics

- **Status**: Wrappers exist for pdpttrsv/pzpttrsv, pddttrsv/pzdttrsv,
  pdpbtrsv/pzpbtrsv, pddbtrsv/pzdbtrsv (interface signatures match
  upstream — UPLO/TRANS/etc.).
- **Symptom**: A naïve "verify pdpttrsv(L) ∘ pdpttrsv(U) ≡ pdpttrs"
  test fails by 60-200x — the trsv routines do NOT compose to the
  full LDL^T solve; the implicit D^-1 step lives only in the *trs
  driver, not in the half-solves.
- **Action**: Build a meaningful test by applying *trsv to a vector
  built from the factor (forward-substitute against the L stored in
  AF), or by rebuilding T from D/E/AF on rank 0 and checking
  T_L * X_got ≡ B_orig directly. Deferred until an upstream-doc
  read pins down exactly what each *trsv variant computes.

## pdposvx / pzposvx — fail via internal pdpocon / pzpocon

- **Status**: pdposvx and pzposvx wrappers exist. Drivers reliably
  abort with "On entry to P{Q,Y}POCON parameter number 10 had an
  illegal value", and on kind16 the gathered X disagrees with the
  reference by 10^3. pdgesvx / pzgesvx pass cleanly on all targets.
- **Diagnosis**: both *posvx drivers funnel through pd/pzpocon for
  rcond. Param 10 is the work pointer; the workspace size that
  *posvx passes into *pocon does not match upstream's expectations
  (or upstream's pocon workspace contract is wrong). Same root
  cause as the existing pdpocon/pzpocon entry above.
- **Action**: Fix or sandbox pd/pzpocon workspace handling, then
  re-enable a *posvx driver. Wrappers remain exposed.

## pdtrsen — heap corruption on every call

- **Status**: target_pdtrsen wrapper exists and the precision result
  is correct (sorted eigenvalue spectrum agrees bit-equally with
  LAPACK dtrsen), but every PDTRSEN call leaves the heap in a
  corrupted state — `free(): invalid pointer` aborts on the next
  deallocate. Same family as pdsyevx (heap corruption via internal
  workspace bookkeeping).
- **Action**: Either fix the upstream PDTRSEN workspace handling
  or run pdtrsen in a sandbox so a single-call result can be reported
  without the cleanup crash. Wrapper remains exposed; pdtrord (the
  reorder-only sibling) does NOT exhibit the issue and ships a
  passing driver.

## pdormrz / pzunmrz update — semantic mismatch with dormrz/zunmrz

- Followup attempt with SIDE='R' (which side-steps the SIDE='L'
  PBETRAN parameter-11 abort) produces a numerically clean run
  (no abort) but the gathered C disagrees with LAPACK dormrz by
  factors of ~1.3 — same magnitude regardless of TRANS='N' vs 'T'.
  Either pdormrz stores the reflectors with a different sign/scale
  convention than dormrz, or the K/L semantics differ. Defer until
  a careful upstream-doc walkthrough.
