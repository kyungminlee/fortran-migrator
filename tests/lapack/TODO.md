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

## Remaining user-facing routines without test drivers (as of session end)

The following user-facing entries from `tests/RESULT.md` still need test drivers.
They are grouped by the original phase plan in `~/.claude/plans/i-want-to-merge-dynamic-crab.md`.

### P8/P9/P10 — Bunch–Kaufman variants + Z-symmetric fillers

D-only:
  sytrf_aa, sytrf_aa_2stage, sytrf_rk, sytrf_rook,
  sytf2_rk, sytf2_rook,
  sytri, sytri2x, sytri_3x, sytri_rook, sycon_rook,
  sytrs_aa, sytrs_aa_2stage, sytrs_rook, syswapr,
  syconv, syconvf, syconvf_rook,
  sysv_aa, sysv_aa_2stage, sysv_rk, sysv_rook,
  sytrd_2stage, sytrd_sy2sb, sb2st_kernels

Z-Hermitian (mirror set):
  hetrf_aa, hetrf_aa_2stage, hetrf_rk, hetrf_rook,
  hetf2_rk, hetf2_rook,
  hetri, hetri2x, hetri_3x, hetri_rook, hecon_rook,
  hetrs_aa, hetrs_aa_2stage, hetrs_rook, heswapr,
  syconv (Z), syconvf (Z), syconvf_rook (Z),
  hesv_aa, hesv_aa_2stage, hesv_rk, hesv_rook,
  hetrd_2stage, hetrd_he2hb, hb2st_kernels

Z-symmetric fillers:
  zsymv, zsyr, zspmv, zspr, zsysv, zsyrfs, zsytrf, zsytri, zsytrs, zsycon,
  zhpcon, zhprfs, zhpsvx, zspsvx

Pattern is mostly a clone-and-modify of the existing `dsytrf` / `dsytrs` /
`dsycon` / `dsyrfs` tests with the appropriate uplo/_aa/_rk/_rook variant
specifier.  The test setup uses a deterministic SPD or Hermitian PD matrix
and compares the ipiv vector + factor element-wise; both implementations
follow the same pivot rule on inputs without ties.

### P17 — Extra-precise expert drivers + gtsv/ptsv family

  d/z × gbsvxx, gesvxx, posvxx, gerfsx, gbrfsx, porfsx, syrfsx, gtsvx, ptsvx, trsna
  d × gtsv, d × ptsv

The xx variants (`gesvxx`, `gerfsxx`, etc.) are the extra-precise iterative
refinement family and are likely to be sensitive on `multifloats`; expect
TODO entries during implementation.

### P22 — Dynamic mode decomposition

  d/z × gedmd, gedmdq

DMD has a wide parameter surface (jobs/jobz/jobr/jobf/whtsvd) -- start with
WHTSVD=1, JOBS='N', JOBZ='V', JOBR='F', JOBF='N' for a basic eigenvalue check.

### Audit pending

After all phases land, run the audit:

```bash
awk 'NR>=89 && NR<=1107 && /^\| ✓ \|/ {gsub(/^\| ✓ \| /,""); split($0,a,/ \|/); print a[1]}' tests/RESULT.md \
  | while read r; do
      ls tests/lapack/*/test_${r}.f90 >/dev/null 2>&1 || echo "MISSING: $r"
    done
```

Expected output: empty (all user-facing routines covered).
