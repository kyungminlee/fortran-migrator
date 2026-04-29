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
