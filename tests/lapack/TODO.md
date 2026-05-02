# tests/lapack — TODO

Issues observed while adding LAPACK tests that need fixes outside
`tests/lapack/` (migrator, build infra, etc.). Each entry names the
specific routines blocked, the file(s) involved, and what would have
to change.

## Phase L5 — *gesvj returns Infinity on multifloats target only — RESOLVED

Status: PASS on all three targets (kind10/kind16/multifloats),
~31 digits on multifloats.

Root cause: upstream multifloats's Blue's-style scaling constants
were defined inverted-in-magnitude: `DD_SBIG = 1e+50` (a SCALE-UP
factor) and `DD_SSML = 1e-50` (a SCALE-DOWN factor) — opposite to
LAPACK's convention where `sbig << 1` (scale-down for big values)
and `ssml >> 1` (scale-up for small values). Compounding that, the
1e±100 / 1e±50 magnitudes only kept `(ax * sbig)**2` representable
for `ax` up to ~1e+104. dgesvj's protect-from-underflow up-scaling
pass routinely produces column norms up to ~`huge(double)/sqrt(N)`
(verified at ~4.5e+153 for the test sizes), so `(ax * sbig)**2`
overflowed to `+Inf` on every column, and the resulting `SVA`
accumulator carried Inf into the Ritz/eigvalue paths.

Fix shape (touches three places):

1. `recipes/lapack/mf_helpers/la_constants_mf.f90` — `mtsml/mtbig/
   mssml/msbig` now compute the LAPACK formula at compile time from
   the leading-limb double's `radix/minexponent/maxexponent/digits`,
   instead of re-exporting upstream multifloats's `DD_*` (which were
   wrong). Fixes the LAPACK side via `USE LA_CONSTANTS_MF` aliases
   in mlassq / wlassq / etc.

2. `recipes/blas/source_overrides/{dnrm2,dznrm2}.f90` (new) —
   replace upstream BLAS's inline `parameter` declarations of
   `tsml/tbig/ssml/sbig` (which the migrator's nuke-rename pass
   substitutes with imports of broken `DD_*`) with renamed locals
   `btsml/btbig/bssml/bsbig` initialized at first call from the
   LAPACK formula (SAVE-vars init'd in an `if (.not. initialized)`
   block — multifloats's `**` operator is user-defined and not
   allowed in PARAMETER initializer expressions). The renamed
   identifiers also dodge the migrator's nuke-rename match.

3. `src/pyengine/fortran_migrator.py` — minor migrator fix in
   `replace_intrinsic_calls` and `replace_generic_conversions`:
   `REAL(x, wp)` was being wrapped to `real64x2(x, wp)` (a 2-arg
   structure constructor — fails because `real64x2` has only one
   component, `limbs`). Now drops the kind-spec second argument
   so the result is `real64x2(x)`, which dispatches through the
   multifloats generic interface (e.g. `dd_from_int`).

| Target      | test_dgesvj      | test_zgesvj      |
|-------------|------------------|------------------|
| kind16      | 2/2 PASS         | 2/2 PASS         |
| kind10      | 2/2 PASS         | 2/2 PASS         |
| multifloats | 2/2 PASS (~31d)  | 2/2 PASS (~31d)  |

Reproduce:
```bash
cd /home/kyungminlee/Code/fortran-migrator/src
uv run python -m pyengine stage /tmp/stg-mf --target multifloats --libraries blas lapack
cmake -S /tmp/stg-mf -B /tmp/stg-mf/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg-mf/build -j8 --target test_dgesvj test_zgesvj
/tmp/stg-mf/build/tests/lapack/test_dgesvj
```

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

## Phase P6/P7 — zgesvdq runtime crash — RESOLVED

Status: PASS on all three targets (kind10/kind16/multifloats),
~31-34 digits.

Root cause: same migrator bug as L5 — `REAL(x, wp)` /
`REAL(x, KIND=wp)` in upstream zgesvdq's body (e.g.
`SQRT(REAL(M, KIND=WP))` for column-norm scaling) was being wrapped
to a 2-arg structure constructor `real64x2(x, wp)` on multifloats,
and the analogous wrap on kind16 sometimes corrupted an `INFO` /
workspace-size argument because the second arg landed in the
wrong slot. Fix landed in `replace_intrinsic_calls` and
`replace_generic_conversions` (drop the kind-spec second argument
when wrapping to the multifloats constructor / when emitting
`REAL(x, KIND=...)` for kind targets). See P22 / L5 entries for
details.

| Target      | test_zgesvdq      |
|-------------|-------------------|
| kind16      | 2/2 PASS (~34d)   |
| kind10      | 2/2 PASS (~18d)   |
| multifloats | 2/2 PASS (~31d)   |

Reproduce:
```bash
cd /home/kyungminlee/Code/fortran-migrator/src
uv run python -m pyengine stage /tmp/stg-q --target kind16 --libraries blas lapack
cmake -S /tmp/stg-q -B /tmp/stg-q/build && cmake --build /tmp/stg-q/build -j8 --target test_zgesvdq
/tmp/stg-q/build/tests/lapack/test_zgesvdq
```

## Phase P22 — complex DMD precision divergence (kind16 / multifloats) — RESOLVED

Status: bit-exact PASS on kind16 / kind10 / multifloats after fixing the
**reference** path. The migrated `xgedmd` / `xgedmdq` were correct all
along; the divergence was a latent fp64 leak in the quad-promoted
reference.

Root cause: `external/lapack-3.12.1/SRC/zgedmd.f90:1006`

```fortran
W(i,j) = CMPLX(RWORK(N+i),ZERO,KIND=WP)*W(i,j)   ! WP = real64 = 8
```

Under `-freal-8-real-16`, gfortran 13 aliases `real(KIND=8)` *storage*
to 16 bytes, but the explicit `CMPLX(...,KIND=8)` intrinsic still
rounds through fp64 internally before storing into the (aliased)
16-byte slot. That bleeds ~1 ULP at fp64 into the Rayleigh quotient,
which the eigenvalue solver then exposes. The migrator rewrites
`WP = real64` to `WP = 16` so the migrated form `CMPLX(...,KIND=16)`
hits the truly-quad branch and is correct. `zgedmdq` inherits the
leak by calling `zgedmd` recursively.

Fix: a reference-only override at `tests/lapack/reflapack/overrides/zgedmd.f90`
that uses `RWORK(N+i)*W(i,j)` (mathematically identical) to avoid the
intrinsic. Wired into the reflapack build via a generic overrides glob
in `tests/lapack/reflapack/CMakeLists.txt` — files in `overrides/`
displace their same-named upstream counterparts in the source list.
Upstream `external/` is untouched.

| Target      | test_dgedmd | test_zgedmd | test_dgedmdq | test_zgedmdq |
|-------------|-------------|-------------|--------------|--------------|
| kind16      | 48/48 PASS  | 48/48 PASS  | 32/32 PASS   | 32/32 PASS   |
| kind10      | 48/48 PASS  | 48/48 PASS  | 32/32 PASS   | 32/32 PASS   |
| multifloats | 48/48 PASS  | 48/48 PASS  | 32/32 PASS   | 32/32 PASS   |

Reproduce:
```bash
cd /home/kyungminlee/Code/fortran-migrator/src
uv run python -m pyengine stage /tmp/stg-q --target kind16 --libraries blas lapack
cmake -S /tmp/stg-q -B /tmp/stg-q/build
cmake --build /tmp/stg-q/build -j8 --target test_zgedmd test_zgedmdq
/tmp/stg-q/build/tests/lapack/test_zgedmd
/tmp/stg-q/build/tests/lapack/test_zgedmdq
```

Test design notes (for future tightening):

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

### P17xx — Extra-precise iterative-refinement family (18 routines) — RESOLVED

  d/z × gbrfsx, gbsvxx, gerfsx, gesvxx, porfsx, posvxx, syrfsx, sysvxx
  z × herfsx, hesvxx

Status: 18/18 PASS on kind10 / kind16 / multifloats (~17–19 digits on
kind10, ~32–34 on kind16, ~30–32 on multifloats). Tests live at
`tests/lapack/linear_solve/test_{d,z}{gesvxx,gerfsx,gbsvxx,gbrfsx,
posvxx,porfsx,sysvxx,syrfsx}.f90` plus `test_zhesvxx.f90` /
`test_zherfsx.f90`. Each driver compares X / RCOND (and RPVGRW for the
SVXX drivers); BERR / ERR_BNDS_NORM / ERR_BNDS_COMP are computed on
both sides but not differenced (they're algorithm-internal estimates,
not deterministic across precision targets).

Resolution had three load-bearing pieces:

1. **Migrated lapack now PUBLIC-links migrated xblas**
   (`cmake/CMakeLists.txt`). The xx-family inside qlapack/elapack/
   mlapack reaches the migrated `qla_*_extended` / `ela_*_extended` /
   etc., which call the migrated XBLAS L2 entry points
   (`BLAS_qgemv_x_`, `BLAS_egemv_x_`, ...). Adding xblas to the lapack
   archive's interface link makes those symbols resolve transitively
   in every test executable that pulls in lapack.

2. **Quad-precision Fortran XBLAS bridge for reflapack_quad**
   (`tests/lapack/reflapack/refxblas_quad_bridge.f90`). The reference
   side `reflapack_quad` is built from the upstream Netlib LAPACK with
   `-freal-8-real-16`, so its `dla_*_extended.f` files call the
   upstream-named `blas_dgemv_x_` / `blas_zhemv_x_` / etc. Those names
   are NOT in qxblas (q-prefixed) and NOT in exblas/mxblas (different
   working precision). The bridge defines them as Fortran subroutines
   compiled with the same `-freal-8-real-16` ABI, forwarding to the
   quad-promoted refblas_quad routines (DGEMV/ZGEMV/DGBMV/ZGBMV/DSYMV/
   ZHEMV) for everything with a BLAS analogue and hand-coding the
   complex-symmetric matvec (no upstream BLAS analog). The Fortran-
   callable signature OMITS the C interface's leading `order` arg
   (the f2c bridge in `external/xblas-1.0.248/src/<routine>/
   BLAS_<routine>_x-f2c.c` hardcodes `blas_colmajor`); getting this
   wrong manifests as `xerbla DGEMV parameter 6 (LDA) had an illegal
   value` because all subsequent args shift by one slot.

3. **Wrapper module entries** in
   `tests/lapack/common/target_lapack_body.fypp` and matching
   interfaces in `tests/lapack/common/ref_quad_lapack.f90`. NPARAMS=0
   in every test driver — disables the PARAMS knob so the upstream
   default (refinement on, threshold 10.0, normwise-only error
   bounds) applies. The 4-call workspace footprint is `4*N` real for
   d-side, `2*N` complex + `2*N` real for z-side, except `zhesvxx`
   which needs `5*N` complex.

Build wiring: `cmake/CMakeLists.txt` makes `${LIB_PREFIX}lapack`
PUBLIC-link `${LIB_PREFIX}xblas` whenever both targets exist. No
opt-in build flag — XBLAS is always built and the bridge is small
(~280 lines of Fortran).

The original blocker description (preserved below) is now obsolete.

#### Original blocker description (for archaeology)

The xx family in LAPACK delegates to `xla_*_extended` helper routines,
which in turn call XBLAS extended-precision BLAS
(`blas_zhemv2_x_`, `blas_zhemv_x_`, `blas_zsymv2_x_`, `blas_zsymv_x_`,
etc.).  Adding tests for any xx routine surfaced the missing XBLAS
symbols at link time and broke **all** other tests by poisoning
`liblapack_test_target.a`.  This is fixed now — see resolution above.

#### Scope (what was vendored)

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
be a future tightening).  All four variants pass bit-exact on all three
targets after the reference-side override at
`tests/lapack/reflapack/overrides/zgedmd.f90` (see resolved P22
section above).

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
