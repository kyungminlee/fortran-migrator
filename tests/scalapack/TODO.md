# tests/scalapack — known upstream / migrator gaps

## Real-MPI exercise via `linux-impi` preset — 2026-04-30 (RESOLVED on kind16)

`cmake --preset=linux-impi` invokes Intel MPI 2021.18 with
`I_MPI_ADJUST_REDUCE=1` and `I_MPI_ADJUST_ALLREDUCE=1` injected via
`MPIEXEC_PREFLAGS` to skip the impi shortpaths that mishandle
BLACS-registered user MPI ops over `REAL(KIND=16)`. ScaLAPACK tests
went from `0% real signal` (singleton MPI) to **170/170 (100%)** on
the kind16 target with a real 2×2 grid.

The 13 ScaLAPACK tests that originally failed under impi were
resolved by 4 distinct fixes — see commits `cfcf023`,
`67c57d2`, and `e469670`:

* **norms (pdlanhs/pzlanhs)** — upstream IAROW double-advance bug;
  patched in source_overrides. See `doc/UPSTREAM_BUGS.md` entry
  "p?lanhs.f IAROW double-advance".
* **equiv (pzgeequ)** — upstream wrong-axis reduction bug; patched
  in source_overrides. See `doc/UPSTREAM_BUGS.md` entry
  "p?geequ.f column-scale reduction wrong axis". (pdgeequ was
  passing only by coincidence; the same fix applies.)
* **QR-fac (pdggrqf/pzggrqf)** — wrapper TAUB used the row axis
  where TAUB is column-distributed; fixed `loccols(descb)` for
  TAUB and bumped the test pre-allocation to `locn_b`.
* **refine (rfs x6) + expert (svx x2)** — the impi preset was
  missing `I_MPI_ADJUST_ALLREDUCE=1`; BLACS *GAMX2D / IGAMX2D
  paths through `MPIR_Reduce_local` crashed in the intranode
  shortpath on user MPI_Op + REAL(KIND=16) payload.

## Banded test sources — workspace / IPIV size fixes (2026-04-30)

13 banded-family tests went from crash/FAIL to PASS on impi after
three test-side fixes:

1. **PDGBTRF / PZGBTRF IPIV under-allocation** — docs claim IPIV size
   `>= DESCA(NB)`, but the divide-and-conquer merge step at
   `pdgbtrf.f:1045` writes `IPIV(LN+1 .. LN+BM+BMN)` with `LN <= NB-BW`
   and `BM+BMN` up to `2*(BW+BWU)`. Tests now allocate
   `nb + 2*(bwl+bwu)`. (Affects `test_p[dz]gbtrf`, `test_p[dz]gbtrs`,
   `test_p[dz]gbsv`.)

2. **PDDBTRS / PZDBTRS workspace under-allocation** — upstream
   `WORK_SIZE_MIN = max(bwl,bwu)*nrhs` is too small; the local-update
   `QLAMOV` at `pqdbtrsv.f:1500` writes `BWU` rows × `NRHS` cols at
   stride `MAX_BW+BWL` starting from `WORK(1+MAX_BW-BWU)`. Tests now
   override `lwork` to `nrhs*max(bwl,bwu) + (nrhs-1)*(bwl+bwu)`.
   (Affects `test_p[dz]dbtrf`, `test_p[dz]dbtrs`, `test_p[dz]dbtrsv`.)

3. **test_pdpttrf invalid reference comparison** — original test
   compared `d_loc, e_loc` against LAPACK `dpttrf`. PDPTTRF's modified
   d/e are NOT LDL^T factors; multi-rank cyclic reduction stores
   reduced-system fillin in `af`, leaving d/e in a representation only
   PDPTTRS understands. Test now does factor+solve and compares solve
   result, mirroring `test_pzpttrf` (which was already correct).

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

## pdgecon / pdpocon (condition-number estimators) — RESOLVED 2026-05-01

Implemented the κ-true reference helper and Hager-Higham bound check.
`tests/scalapack/common/cond_helpers.f90` exposes
`true_kappa1_general` / `true_kappa1_general_z` (LU + dgetri / zgetri
inverse) and `true_kappa1_posdef` / `true_kappa1_posdef_z` (Cholesky +
dpotri / zpotri). `tests/scalapack/auxiliary/test_pdgecon.f90` and
`test_pdpocon.f90` factor the matrix on the distributed side, call the
estimator, and assert
``rel_err_scalar(kappa_est, kappa_true) <= 0.7`` — covers the
worst-case 1/3 LACON ratio with margin for fp noise.

Results (kind16, 2×2 grid):
- `test_pdgecon`  3/3 PASS, κ_est ≈ 0.78 · κ_true (~22% rel-err)
- `test_pdpocon`  3/3 PASS, κ_est ≈ 0.83 · κ_true (~17% rel-err)

(pzheev fix folded into the pxhetrd/pxheev entry above.)

## pdormrz / pzunmrz — PARTIALLY RESOLVED 2026-05-02 (SIDE='R' only)

Two upstream bugs identified and patched in
`recipes/scalapack/source_overrides/p[dz]larzb.f` and
`recipes/scalapack/source_overrides/p[dz]ormrz.f`:

1. **`p[dz]larzb.f` PBDTRAN N-arg vs LV-buffer mismatch (SIDE='L'
   branch).** With `mb=nb=8`, `IA=JA=IC=JC=1`, `K=mA`, `L=nA-mA`, and
   `mC=nA=K+L` (the only legal shape per `pdormrz.f:34` "Q is of order
   M if SIDE='L'"), `pdlarzb.f:374-377` allocates `WORK(IPV)` with
   `LV = MAX(1, MPC20) = local rows of L-row matrix`, but
   `PBDTRAN(...,'Rowwise','Transpose', K, M+ICOFFV, ..., LV, ...)`
   checks `LDC ≥ NP` where `NP = local rows of M=mC=(K+L)-row matrix`.
   Since `mC > L` for any non-trivial K, `LV < NP` and `PBDTRAN` exits
   via `PXERBLA(INFO=11)` without performing the transpose, leaving
   uninitialised `WORK(IPV)`. The override changes the PBDTRAN N arg
   from `M+ICOFFV` to `L+ICOFFV` — only the trailing K×L portion of V
   holds meaningful Householder data per the row-stored RZ convention.

2. **`p[dz]ormrz.f` post-loop condition copy-paste (`(L&&!T)||(!L&&T)`
   matches the pre-loop condition instead of being its complement).**
   The result is that for SIDE/TRANS combinations `(LEFT,NOTRAN)` and
   `(RIGHT,TRANS)`, neither pre nor post fires. Combined with the main
   DO loop's asymmetric I1/I2 bounds, the leading partial block of
   reflectors goes unapplied. The override fixes the post-loop
   condition to `(LEFT && NOTRAN) || (.NOT.LEFT && .NOT.NOTRAN)`,
   mirroring `pdormrq.f:454`'s correct shape.

Both fixes wired in `recipes/scalapack.yaml` (`source_overrides` +
`prefer_source` pins).

**SIDE='R' status**: passes to target precision on kind16 / 2×2 grid
for both TRANS='N' and TRANS='C/T'. Test drivers:
`tests/scalapack/factorization/test_p[dz]ormrz.f90`.

**SIDE='L' status — STILL DEFERRED**. Even with both fixes applied,
SIDE='L' (both TRANS='N' and TRANS='T') fails by a factor of ~2.5.
The failure reproduces with K=mA=4 / mC=nA=K+L=8 (single-block path,
where PDLARZB is *not* called and only the post-loop PDORMR3 fires),
ruling out PDLARZB. The bug appears to live in the PDORMR3 + PDLARZ
chain for SIDE='L' (the SIDE='R' chain works correctly). Test drivers
currently restrict to SIDE='R' only; re-add SIDE='L' once the
underlying PDORMR3/PDLARZ bug is identified.

Documented in `doc/UPSTREAM_BUGS.md`.

## *trsv (single-half banded/tridiag triangular solves) — RESOLVED 2026-05-02

The 8 smoke tests at `tests/scalapack/linear_solve/test_p[dz]{pt,dt,pb,db}trsv.f90`
were upgraded from INFO=0 checks to validation tests that compose the
two halves and compare against the full `*trs` oracle on the same
distributed factor. Per-family composition:

* **PT** (real LDL^T): `pdpttrsv('L')` → mid-step → `pdpttrsv('U')`
  ≡ `pdpttrs`. The mid-step mirrors `pdpttrs.f:707-716` exactly: ranks
  0..NPCOL-2 divide entries 1..nb-1 by `d_loc` (the modified LDL^T
  diagonal D'), then divide entry nb by the reduced-system pivot
  stored at `af(nb+1)`; the last rank divides all nb entries by
  `d_loc` (no boundary). The original "fails by 60-200x" smoke test
  was missing both the D'^-1 step and the AF reduced-system pivot.
* **PT** (complex Hermitian LDL^H): `pzpttrsv('L','N')` →
  mid-step (real D' divide + complex AF pivot divide, mirroring
  `pzpttrs.f:744-758`) → `pzpttrsv('L','C')` ≡ `pzpttrs(UPLO='L')`.
  Note backward half uses `('L','C')` (apply L^H of stored L), not
  `('U','N')` — pzpttrf stores the L view, and `pzpttrs(UPLO='L')`
  internally uses the same pattern.
* **DT/DB** (LU): `trsv('L','N')` ∘ `trsv('U','N')` ≡ `*trs` (no
  scaling). For DB, the `*trsv` and `*trs` paths must use distinct AF
  copies — even though AF is documented read-only, in practice the
  `*trs` divide-and-conquer path perturbs internal scratch in AF on
  certain ranks, so the test snapshots AF after the factor and runs
  the two paths sequentially with restore in between. Workspace must
  be the bumped formula `nrhs*max(bwl,bwu) + (nrhs-1)*(bwl+bwu)`
  (per the banded TODO entry above) — upstream's `*trsv` LWMIN query
  under-allocates.
* **PB** (Cholesky, UPLO='U'): `trsv('U','T'/'C')` then `trsv('U','N')`
  ≡ `pbtrs` — apply U^-T (or U^-H) before U^-1.

All 8 PASS to ~target precision on kind16 / 2×2 grid.

## pdposvx / pzposvx — RESOLVED 2026-05-01

Two upstream LWMIN bugs in `pdposvx.f` / `pzposvx.f`. Both were
fixed via `recipes/scalapack/source_overrides/p[dz]posvx.f` plus
`prefer_source: PDPOSVX, PZPOSVX` in `recipes/scalapack.yaml`.

1. `LWMIN = 3*DESCA(LLD_)` is too small for the internal PDPOCON /
   PZPOCON contract (param-10 abort). Documented contract is
   `MAX(PDPOCON_LWMIN, PDPORFS_LWMIN)`.
2. Companion bug in `pzposvx.f`: `LRWMIN = MAX(2*NQ, NP)` doesn't
   cover PZLANHE('1', ...) which needs `2*NQMOD + NPMOD + LDW`.
   PZLANHE then writes past RWORK and corrupts the heap, surfacing
   as a `free()` SIGSEGV during cleanup.

Both bugs documented in `doc/UPSTREAM_BUGS.md`. Tests:
`tests/scalapack/linear_solve/test_pdposvx.f90` / `test_pzposvx.f90`
both 3/3 PASS at ~33-digit accuracy on kind16 / 2×2 grid.

## pdtrsen — RESOLVED 2026-05-02

Root cause was an upstream LQUERY-contract bug: `pdtrsen.f:499-538`
writes `IWORK(1:N)` (the SELECT→integer copy plus the IGAMX2D
broadcast) before the LQUERY return at line 619-622, so a caller
passing `IWORK(1)` for the workspace query writes past the end of its
buffer and corrupts the heap. The follow-up "real" call runs cleanly
(LIWMIN ≥ N), which is why the eigenvalue spectrum was correct even
as the heap was already damaged.

Fix: `target_pdtrsen` in
`tests/scalapack/common/target_scalapack_body.fypp` allocates a local
`iwork_t(max(1,n))` in the LQUERY branch and forwards it to upstream
in place of the caller's IWORK; `iwork_t(1)` (= LIWMIN) is copied back
to `iwork(1)` afterwards. Mirrors the `WORK_T(3)` pattern for
`target_pdsyevx`'s upstream `WORK(1:3)` early-write.

New driver `tests/scalapack/eigenvalue/test_pdtrsen.f90` exercises the
fix over sizes [32, 64, 96]; all PASS without the previous `free()`
abort. Documented in `doc/UPSTREAM_BUGS.md`.

## pdormrz / pzunmrz "factor of ~1.3" — folded into the resolution above (2026-05-02)

Same root cause as the SIDE='L' bug: `pdlarzb.f`'s PBDTRAN N-arg vs
LV-buffer mismatch. The transpose silently no-ops via PXERBLA, leaving
uninitialised workspace that gets multiplied through to C —
manifesting as a ~1.3× scaling error. Resolved by the same
`source_overrides` fix.
