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

## Phase P22 — complex DMD precision divergence (kind16 / multifloats)

`tests/lapack/eigenvalue/test_zgedmd.f90` and
`tests/lapack/eigenvalue/test_zgedmdq.f90` fail layer-(a) (sorted Ritz
spectrum) on kind16 and multifloats: the migrated target (`xgedmd` /
`xgedmdq`) disagrees with the quad reference (`zgedmd` / `zgedmdq`
compiled with `-freal-8-real-16`) by ~1-5 ULPs at *double* precision
(max_rel_err 1.9e-17 .. 4.6e-17). On kind10 the tolerance window absorbs
the divergence and tests pass. Layer-(b) (residual differential) tracks
the layer-(a) divergence — fails on the same targets, passes on kind10.

| Target      | test_dgedmd | test_zgedmd | test_dgedmdq | test_zgedmdq |
|-------------|-------------|-------------|--------------|--------------|
| kind16      | 48/48 PASS  | 12/48       | 32/32 PASS   | 0/32         |
| kind10      | 48/48 PASS  | 48/48 PASS  | 32/32 PASS   | 32/32 PASS   |
| multifloats | 48/48 PASS  | 12/48       | 32/32 PASS   | 0/32         |

Real siblings `test_dgedmd` and `test_dgedmdq` pass cleanly on all
three targets, so the bug is specific to the complex prefix migration.

Hypothesis: a fp64-precision step inside `xgedmd` / `xgedmdq` not present
in `xgesvd` / `xgeev` themselves (those are bit-exact on kind16 against
their references — verified via `test_zgesvd` / `test_zgeev`). The
migration diff between `external/lapack-3.12.1/SRC/zgedmd.f90` and
`lapack/src/xgedmd.f90` is mechanical (Z→X / D→Q / DBLE→REAL(.,KIND=16))
with no obvious fp64 leak — the bug likely lives in one of the migrated
helpers transitively reachable from `xgedmd` / `xgedmdq` but not from
the leaf SVD/eigenvalue routines (candidates: `xlassq` chain, `xlascl`
for the Rayleigh-quotient build, or a CMPLX/REAL-style intrinsic call
that resolves to the wrong kind in the complex path).

Reproduce:
```bash
cd /home/kyungminlee/Code/fortran-migrator/src
uv run python -m pyengine stage /tmp/stg-q --target kind16 --libraries blas lapack
cmake -S /tmp/stg-q -B /tmp/stg-q/build
cmake --build /tmp/stg-q/build -j8 --target test_zgedmd test_zgedmdq
/tmp/stg-q/build/tests/lapack/test_zgedmd
/tmp/stg-q/build/tests/lapack/test_zgedmdq
```

Tests left in place per the failing-test-stays-visible convention.

Test design notes (for future tightening once the bug is fixed):

  - Sizes for `*gedmdq` capped at `(12, 8)` (two-size sweep instead of
    the design's `(24, 16)`): single-trajectory snapshot conditioning
    `κ(F) ∝ |λ_max/λ_min|^N` collapses fast, so larger N loses
    precision on kind10/multifloats even with relaxed tolerance.
  - Layer-(a) tolerance: `100 * n^2 * eps` for `*gedmd`, `10000 * n^2 * eps`
    for `*gedmdq` (trajectory data needs the bump).
  - Layer-(b) is the differential `‖res_got − res_ref‖∞ / max(‖A‖_F, eps)`
    rather than a relative-to-zero ratio, which collapses on near-zero
    residuals (multifloats target manifests this most clearly).
  - Layer-(c) (reconstruction `‖Y − Z·diag(λ)·Z⁺·X‖_F / ‖Y‖_F`) is not
    implemented in v1; it would catch additional Z-discrepancy modes.

## Remaining user-facing routines without test drivers (as of 2026-04-29)

The Bunch–Kaufman / Z-symmetric / packed / utility / 2-stage tridiag (incl.
sy2sb/he2hb) / gtsvx / ptsvx / trsna / DMD (gedmd/gedmdq) families are now
covered (Phases P8a–P8h, P22).  20 routines remain — they cluster into
two categories, both of which are blocked or deferred:

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

### P22 — Dynamic mode decomposition (4 routines) — DELIVERED

  d/z × gedmd, gedmdq

Tests landed in `tests/lapack/eigenvalue/test_{d,z}gedmd{,q}.f90` plus
helpers (`gen_orthogonal_quad` in `test_data.f90`, `sort_eig_lex_z` in
`compare.f90`), 4 ref interfaces in `ref_quad_lapack.f90`, and 4 fypp
wrappers in `target_lapack_body.fypp`.  Layer-(a) spectrum + layer-(b)
differential residuals; layer-(c) reconstruction not implemented (would
be a future tightening).  Real variants pass on all three targets;
complex variants surface the migrator-side fp64 leak documented in the
"Phase P22 — complex DMD precision divergence" section above.

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
