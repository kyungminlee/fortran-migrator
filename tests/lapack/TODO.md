# tests/lapack — TODO

Issues observed while adding LAPACK tests that need fixes outside
`tests/lapack/` (migrator, build infra, etc.). Each entry names the
specific routines blocked, the file(s) involved, and what would have
to change.

## Phase L5 — *gesvj returns Infinity on multifloats target only

`tests/lapack/svd/test_dgesvj.f90` and `test_zgesvj.f90` pass on
`kind16` and `kind10` but fail on `multifloats` — the target's
`tgesvj`/`vgesvj` returns `Infinity` for all singular values while
the quad reference returns sensible values. Other Jacobi-SVD routines
(`dgejsv`/`zgejsv`) pass on all targets.

Hypothesis: the migrated `t/vgesvj` may mishandle the default
convergence-tolerance path (`ctol=0` → routine selects own threshold)
in the `real64x2` arithmetic — possibly in a `LAMCH`-like inquiry or
in an internal `MAX/MIN` call that returns NaN/Infinity for the dd type.

To reproduce:
```bash
cd /home/kyungminlee/Code/fortran-migrator/src
uv run python -m pyengine stage /tmp/stg-mf --target multifloats --libraries blas lapack
cmake -S /tmp/stg-mf -B /tmp/stg-mf/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg-mf/build -j8 --target test_dgesvj test_zgesvj
/tmp/stg-mf/build/tests/lapack/test_dgesvj
```

Tests left in place per the failing-test-stays-visible convention.

## Phases L6..L23 — not yet started

Remaining phases per `~/.claude/plans/start-a-project-to-stateless-bumblebee.md`:

  L5  Modern least-squares + Jacobi SVD (gelsd/gelss/gelsy/gelst/gejsv/gesvj)
  L6+L7 Modern SVD + GSVD + CSD (gesvdq/gesvdx/bdsvdx/bbcsd/ggsvd3/ggsvp3/orcsd/orcsd2by1)
  L8  Bunch-Kaufman variants — factorization (sytrf_rk/_rook/_aa/_aa_2stage/sytri_3/_rook/syconv*/syconvf*)
  L9  Bunch-Kaufman variants — solve/inverse/cond (sytri/sytri2/sytri2x/sytrs2/sytrs_3/_aa/_aa_2stage/_rook/sycon_3/_rook)
  L10 Bunch-Kaufman driver families (sysv_rk/_rook/_aa/_aa_2stage/sysvxx/syswapr)
  L11 Pivoted Cholesky + RFP (pstrf/pbstf/pteqr/pftrf/pftri/pftrs/tftri/sfrk)
  L12 Storage conversion RFP/tri/packed (tfttr/tfttp/tpttr/tpttf/trttp/trttf) — ship before L11
  L13 Pentagonal QR/LQ (tpqrt/_2/tpmqrt/tplqt/_2/tpmlqt/tprfb)
  L14 Householder reconstruction & TSQR (orhr_col/orgtsqr/_row/getsqrhrt/geqp3rk/tzrzf)
  L15 CS decomposition (orbdb/_1/_2/_3/_4/_5/_6/orm22)
  L16 Tridiagonal eigensolvers MRRR (stegr/stein/stemr)
  L17 Expert drivers extra-precise XX (gesvxx/gbsvxx/posvxx/gerfsx/gbrfsx/porfsx/syrfsx/gtsvx)
  L18 Auxiliary norm utilities (langb/langt/lanhs/lansb/lansf/lansp/lanst/lansy/lantb/lantp/lantr)
  L19 Permutation/norm helpers (lapmr/lapmt/lapll/lacn2/lacon/lartg/lartgp/lartgs)
  L20 Small public scalar utilities (lapy2/lapy3/ladiv/lamrg/larnv)
  L21 Generalized symmetric/Hermitian glue (sbgst/spgst/sygst/sbgvx/spgvx/sygvx)
  L22 Modern dynamic mode decomposition (gedmd/gedmdq)
  L23 Audit / consolidation

## Older blocked-T routines that need a wider scope to test

The following routines were tested in Phase 35 / Phase 36, but only
the |R| (or |L|) factor is compared — the lower-triangle reflector
storage and the workspace `T` (whether the reserved workspace of
`dgeqr` or the explicit block-T of `dgeqrt`) differ between
implementations because the recursion / blocking choices in the
migrated and reference paths diverge:

- `dgeqr`/`zgeqr`, `dgelq`/`zgelq` (Phase 35)
- `dgeqrt`/`zgeqrt`, `dgelqt`/`zgelqt` (Phase 36)

For the matching `*gemqr`/`*gemlq` and `*gemqrt`/`*gemlqt` apply
routines the C output is canonical and matches cleanly — those tests
do compare the full result.

If migrator changes (or recipe overrides) ever pin both paths to the
same blocking heuristic, the factor-side comparisons could be
strengthened to compare the full A and T arrays.

## Phase P6/P7 — zgesvdq runtime crash

`tests/lapack/svd/test_zgesvdq.f90` builds but crashes at runtime on the
kind16 differential test with a Fortran runtime "Expected INTEGER for item 3
in formatted transfer, got CHARACTER" error originating after the call to
`zgesvdq` returns. Suspect: the migrated `xgesvdq` corrupts an integer workspace
size or `info` argument, which then trips a downstream WRITE statement.

Reproduce:
```bash
cd /home/kyungminlee/Code/fortran-migrator/src
uv run python -m pyengine stage /tmp/stg-q --target kind16 --libraries blas lapack
cmake -S /tmp/stg-q -B /tmp/stg-q/build && cmake --build /tmp/stg-q/build -j8 --target test_zgesvdq
/tmp/stg-q/build/tests/lapack/test_zgesvdq
```

dgesvdq passes on the same target; the bug is specific to the complex variant.
Test left in place; needs migrator-side investigation of the
`xgesvdq` call sequence vs `dgesvdq`.

## Remaining user-facing routines without test drivers (as of 2026-04-29)

The Bunch–Kaufman / Z-symmetric / packed / utility / 2-stage tridiag (incl.
sy2sb/he2hb) / gtsvx / ptsvx / trsna families are now covered (Phases
P8a–P8h).  24 routines remain — they cluster into three categories of
non-trivial work:

### P17xx — Extra-precise iterative-refinement family (18 routines) — BLOCKED on XBLAS

  d/z × gbrfsx, gbsvxx, gerfsx, gesvxx, porfsx, posvxx, syrfsx, sysvxx
  z × herfsx, hesvxx

**Blocker:** the xx family in LAPACK delegates to `xla_*_extended` helper
routines, which in turn call XBLAS extended-precision BLAS
(`blas_zhemv2_x_`, `blas_zhemv_x_`, `blas_zsymv2_x_`, `blas_zsymv_x_`,
etc.).  XBLAS is a separate library that the migrator/recipe does not
currently provide for either the `qlapack` (kind16) reference build or
the migrated targets.  Adding tests for any xx routine surfaces the
missing XBLAS symbols at link time and breaks **all** other tests by
poisoning `liblapack_test_target.a`.

Reproduce link error:
```bash
cd /home/kyungminlee/Code/fortran-migrator/src
uv run python -m pyengine stage /tmp/stg-q --target kind16 --libraries blas lapack
# Add even a stub test target_dgesvxx wrapper in target_lapack_body.fypp,
# rebuild liblapack_test_target.a — then any test executable fails to link:
#   undefined reference to `blas_zhemv2_x_'
#   undefined reference to `blas_zsymv_x_'
```

Until XBLAS lands, **do not** add interfaces or wrappers for any xx routine —
the wrapper module pulls in the xla_* objects at link time even if no test
uses the wrapper, breaking the whole test suite.

**Resolution: vendor XBLAS as a separate library, opt-in via a build flag.**

#### Scope (what to vendor)

The xx family has 28 LAPACK helper files that reach into XBLAS:

```
{s,d,c,z}la_{gerfsx,gbrfsx,porfsx,syrfsx}_extended.f      (16 files)
{c,z}la_herfsx_extended.f                                 (2 files)
{s,d,c,z}la_{gercond,gbrcond,porcond,syrcond}_x.f         (8 files — internal only, no XBLAS)
{c,z}la_hercond_x.f                                       (2 files — internal only, no XBLAS)
```

The 16 `*la_*rfsx_extended.f` files plus the 2 Hermitian variants call
exactly **28 XBLAS entry points** — every one is a Level-2 matvec:

```
BLAS_{S,D,C,Z}{GEMV,GBMV,SYMV}_X  + matching 2_X  variants   (24 routines)
BLAS_{C,Z}HEMV_X                  + matching 2_X  variants   (4 routines)
```

(The `*cond_x.f` files only call internal LAPACK helpers — `*la_lin_berr`,
`*la_geamv`, etc. — no XBLAS.)

#### Precision strategy per target

XBLAS does residual refinement in *higher* precision than the working
precision via a (head, tail) accumulator pair.  The migrator promotes the
internal extension precision per target:

| Target       | Working precision | Extension (head, tail)                |
|--------------|-------------------|---------------------------------------|
| kind16       | `real128` (quad)  | `(real128, real128)` — DD-on-quad     |
| kind10       | `real(kind=10)`   | `(real(kind=10), real(kind=10))`      |
| multifloats  | `real64x2` (DD)   | `(real64x2, real64)` — keep DD as     |
|              |                   | head, one extra `real64` for the tail |

Note for multifloats: do **not** introduce a `real64x3` triple-double type.
The (head, tail) decomposition uses the existing DD type for `head` and one
plain `real64` scalar for `tail`.  This matches XBLAS's existing internal
pattern and reuses the DD building blocks already in `multifloats.hh`.

#### Recipe shape

- New `external/xblas-1.0.248/` (or current Netlib XBLAS release) — the
  ~28 entry-point Fortran sources plus the internal helpers they call.
  Public domain license.
- New `recipes/xblas.yaml` mirroring the structure of `recipes/blas.yaml`.
  Goes through the same precision-substitution pipeline as BLAS so each
  target gets its own `libxblas-*.a` (and `libqxblas-*.a` for the kind16
  reference).
- Optional: a sweep of which Netlib XBLAS sources are transitively
  reachable from the 28 entry points; this trims the recipe to the
  minimal vendored subset.

#### Opt-in build flag

Two layers need the same gate:

1. **Recipe / build flag** — e.g. `-DBUILD_XBLAS=ON` at top-level cmake
   configure.  When OFF: skip building libxblas, exclude `*la_*_extended.f`
   from the LAPACK builds (else qlapack still has unresolved XBLAS refs
   even though no test uses them — `liblapack_test_target.a` poisoning
   recurs).  When ON: ship libxblas and include the xla_* objects in
   `liblapack`/`libqlapack`.
2. **Test glob in `tests/lapack/CMakeLists.txt`** — gate the xx-family
   `test_*.f90` files on the same flag, *and* gate the wrapper definitions
   in `target_lapack_body.fypp` (e.g. `#:if BUILD_XBLAS` blocks) so the xx
   `target_*` symbols don't even appear in `liblapack_test_target.a` when
   the flag is off.

Default state of the flag, naming convention, and whether it lives at
top-level or per-recipe are open and should follow whatever the migrator
team prefers.

Once that lands, the xx tests follow the standard pattern (interfaces in
`ref_quad_lapack.f90` + wrappers in `target_lapack_body.fypp` + per-routine
`test_*.f90` files comparing X / RCOND / berr / err_bnds_norm /
err_bnds_comp).  The signatures and partial wrapper code from this
session's aborted attempt are recoverable from the git reflog if useful.

### P22 — Dynamic mode decomposition (4 routines) — pending, rich coverage required

  d/z × gedmd, gedmdq

User has requested **rich coverage** (not minimal smoke).  DMD does not
depend on XBLAS — independent of the P17xx blocker.

#### Routine signatures (LAPACK 3.12.1)

```
DGEDMD  ( JOBS, JOBZ, JOBR, JOBF, WHTSVD,
          M, N, X, LDX, Y, LDY, NRNK, TOL,
          K, REIG, IMEIG, Z, LDZ, RES,
          B, LDB, W, LDW, S, LDS,
          WORK, LWORK, IWORK, LIWORK, INFO )            -- 31 args

DGEDMDQ ( JOBS, JOBZ, JOBR, JOBQ, JOBT, JOBF, WHTSVD,
          M, N, F, LDF, X, LDX, Y, LDY, NRNK, TOL,
          K, REIG, IMEIG, Z, LDZ, RES,
          B, LDB, V, LDV, S, LDS,
          WORK, LWORK, IWORK, LIWORK, INFO )            -- 35 args

ZGEDMD / ZGEDMDQ : same arg lists; eigenvalues come out in EIGS(*) (complex)
                   instead of the (REIG, IMEIG) real pair.
```

#### Sources of harmless ref/target divergence

1. **WHTSVD selection** rotates SVD vectors differently across LAPACK's
   four SVD backends (1=GESVD, 2=GESDD, 3=GESVDQ, 4=GEJSV).  Pin
   WHTSVD=1 across ref and target so both call GESVD.
2. **Eigenvalue ordering** — internal Schur problem returns spectrum in
   implementation-dependent order; compare on a sorted set.
3. **Conjugate-pair sign ambiguity** for the D variant — REIG/IMEIG
   pair order within a complex-conjugate cluster is not canonical.
4. **Eigenvector phase / scaling** — Z columns differ by a unit
   complex scalar (D: real ±1 modulo conjugate-pair convention).

#### Test design (rich coverage)

Three comparison layers per case:

  **(a)** Spectrum check — assemble eigenvalues as `complex` (D: from
  REIG+i·IMEIG; Z: from EIGS), sort lexicographically by (Re, Im) on
  both sides, and compare via `max_rel_err_vec_z`.

  **(b)** Residual norm (RES output) — RES(j) is the LAPACK-reported
  ‖A·z_j − λ_j·z_j‖₂.  Just verify `RES(j) <= tol_res * ‖A‖_F` for
  the target — phase-invariant, no sorting needed.

  **(c)** Reconstruction check (when JOBZ='V', JOBF='N') —
  ‖Y − Z · diag(λ) · (Z⁺ X)‖_F / ‖Y‖_F bounded.  Also phase-invariant.

#### Parameter sweep

For each of `gedmd` and `gedmdq`, run the cross-product:

  - sizes:    `(M, N) ∈ {(8,5), (12,8), (24,16)}`  (m≥n required)
  - JOBS:    `{'N', 'S'}`              (no scaling, column scaling)
  - JOBZ:    `{'N', 'V'}`              (skip Z, full Z)
  - JOBR:    `{'F', 'S'}`              (no residual, residuals)
  - JOBF:    `{'N', 'E'}`              (no exact-DMD, exact-DMD modes)
  - WHTSVD:  `1` only (pinned for cross-impl agreement)
  - NRNK:    `-1` (auto rank from TOL); TOL = 10·target_eps
  - For `gedmdq`: JOBQ ∈ {'F'}, JOBT ∈ {'F'}, JOBF as above

For each combo, apply layers (a)+(b)+(c) where the chosen JOB flags
make the relevant outputs available.  Skip (c) when JOBZ='N'.

#### Tolerances

  - Spectrum (layer a):  `tol_eig = 100 * n^2 * target_eps`
  - Residual (layer b):  `tol_res = 100 * n   * target_eps`
  - Reconstruction (c):  `tol_rec = 100 * n^2 * target_eps`

(Conservative; tighten after first kind16 pass if `err << tol`.)

#### Test data generation

Build A with controlled spectrum:
  `A = U · diag(λ) · U⁻¹` with random orthogonal U from
  `gen_orthogonal_quad`, λ chosen with at least one complex pair
  (D variant) and one near-zero magnitude.
Snapshots: `X = random m×n`, `Y = A · X`.  Ensures the spectrum
the routine recovers is known up to sort order, which is useful for
debugging when (a) fails.

#### File layout

  - `tests/lapack/eigenvalue/test_dgedmd.f90`
  - `tests/lapack/eigenvalue/test_zgedmd.f90`
  - `tests/lapack/eigenvalue/test_dgedmdq.f90`
  - `tests/lapack/eigenvalue/test_zgedmdq.f90`
  - `ref_quad_lapack.f90` — 4 new explicit interfaces
  - `target_lapack_body.fypp` — 4 new `target_*` wrappers (q2t/t2q on
    X, Y, F, output Z, eigenvalues, RES, B, W, V, S as appropriate)

#### Risks / open items

  - `gedmdq` outputs an extra `T` matrix (N×N small upper-triangular
    snapshot) under JOBT='F'; verify whether T is canonical or
    impl-defined.  If non-canonical, skip layer-(c) check that would
    rely on T.
  - `JOBF='E'` returns "exact DMD" Koopman modes in B; sign/phase
    convention may diverge across SVD backends even at WHTSVD=1.
    First implementation should report B with phase-invariant norms,
    not elementwise.
  - Workspace queries: WORK(1)/LWORK and IWORK(1)/LIWORK must be sized
    via the documented LWORK=-1 query path; the wrapper allocates from
    the query result, never a hardcoded floor.

### sb2st kernels (2 routines)

  dsb2st_kernels, zhb2st_kernels — SBR (Successive Band Reduction) inner
  kernels, normally driven by `*sb2st` / `*hb2st`.  Testing them in
  isolation requires constructing a partial band-reduction state that
  matches the kernel's mid-stream invariants — non-trivial, may not be
  worth standalone coverage if the orchestrator paths are tested.

### Audit pending

After all phases land, run the audit:

```bash
awk 'NR>=89 && NR<=1107 && /^\| ✓ \|/ {gsub(/^\| ✓ \| /,""); split($0,a,/ \|/); print a[1]}' tests/RESULT.md \
  | while read r; do
      ls tests/lapack/*/test_${r}.f90 >/dev/null 2>&1 || echo "MISSING: $r"
    done
```

Expected output: empty (all user-facing routines covered).
