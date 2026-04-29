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

Resolution path: extend the recipe (or add an override) to either
(a) vendor XBLAS into `qlapack` and the migrated targets, or
(b) exclude `xla_*_extended.f` from the qlapack/migrated builds and stub
out the xx callers.

Until then, **do not** add interfaces or wrappers for any xx routine —
the wrapper module pulls in the xla_* objects at link time even if no
test uses the wrapper, breaking the whole test suite.

### P22 — Dynamic mode decomposition (4 routines)

  d/z × gedmd, gedmdq

DMD has a wide parameter surface (jobs/jobz/jobr/jobf/whtsvd).  Suggested
starting point: WHTSVD=1, JOBS='N', JOBZ='V', JOBR='F', JOBF='N' for a
basic eigenvalue check.  Not attempted yet; should be doable once the
above XBLAS issue is sorted (DMD does not depend on XBLAS).

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
