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

The Bunch–Kaufman / Z-symmetric / packed / utility / 2-stage tridiag /
gtsvx / ptsvx / trsna families are now covered (Phases P8a–P8f + 2-stage
tridiag).  26 routines remain — they cluster into three categories of
non-trivial work:

### P17xx — Extra-precise iterative-refinement family (18 routines)

  d/z × gbrfsx, gbsvxx, gerfsx, gesvxx, porfsx, posvxx, syrfsx, sysvxx
  z × herfsx, hesvxx

These have ~25-argument signatures including ERR_BNDS_NORM, ERR_BNDS_COMP,
N_ERR_BNDS, NPARAMS, PARAMS for the extra-precise residual refinement.  The
xx routines are also known sensitive on `multifloats`.  Expect to need
careful workspace shapes and per-target tolerance bumps.

### P22 — Dynamic mode decomposition (4 routines)

  d/z × gedmd, gedmdq

DMD has a wide parameter surface (jobs/jobz/jobr/jobf/whtsvd).  Suggested
starting point: WHTSVD=1, JOBS='N', JOBZ='V', JOBR='F', JOBF='N' for a
basic eigenvalue check.

### 2-stage band reductions (4 routines)

  dsytrd_sy2sb, zhetrd_he2hb     — symmetric/Hermitian to banded
  dsb2st_kernels, zhb2st_kernels — band-to-tridiag SBR kernels

`dsytrd_sy2sb` workspace query returns invalid LWORK in the quad reference
build (`ILAENV2STAGE` returns 0 / negative for the kind16 reference
configuration), causing a `Fatal glibc error: malloc.c assertion failed`
when the routine writes out-of-bounds in the workspace setup.  The same
crash reproduces with a hand-sized large workspace (`4*kd*kd*n+4*n*n`),
pointing to a deeper interaction between the migrator output and
`ILAENV2STAGE`/blocking choices for these routines.

`*sb2st_kernels` are the SBR (Successive Band Reduction) inner kernels and
are normally driven by `*sb2st` / `*hb2st`.  Testing them in isolation
requires constructing a partial band-reduction state that matches the
kernel's mid-stream invariants — non-trivial, may not be worth standalone
coverage if the orchestrator paths are tested.

Reproduce sy2sb crash:
```bash
cd /home/kyungminlee/Code/fortran-migrator/src
uv run python -m pyengine stage /tmp/stg-q --target kind16 --libraries blas lapack
cmake --build /tmp/stg-q/build -j8 --target test_dsytrd_sy2sb   # would require restoring removed test
```

### Audit pending

After all phases land, run the audit:

```bash
awk 'NR>=89 && NR<=1107 && /^\| ✓ \|/ {gsub(/^\| ✓ \| /,""); split($0,a,/ \|/); print a[1]}' tests/RESULT.md \
  | while read r; do
      ls tests/lapack/*/test_${r}.f90 >/dev/null 2>&1 || echo "MISSING: $r"
    done
```

Expected output: empty (all user-facing routines covered).
