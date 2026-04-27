# tests/lapack — TODO

Issues observed while adding LAPACK tests that need fixes outside
`tests/lapack/` (migrator, build infra, etc.). Each entry names the
specific routines blocked, the file(s) involved, and what would have
to change.

## 2-stage eigenvalue routines: untestable as currently migrated

**Blocked tests:** `dsyev_2stage`, `zheev_2stage`, `dsyevd_2stage`,
`zheevd_2stage`, plus the related `dsbev_2stage`, `dsbevd_2stage`,
`dsbevx_2stage`, `dsyevr_2stage`, `dsyevx_2stage`, and complex
counterparts.

Three compounding gaps:

1. **`.F` files excluded from manifest globs.**
   `src/pyengine/__main__.py` (around line 770) globs
   `*.f`, `*.f90`, `*.F90` but not `*.F`. The 2-stage chain depends on
   `iparam2stage.F`, `dsytrd_sb2st.F`, `zhetrd_hb2st.F` (capital-F,
   OpenMP-guarded). Without `*.F` they never enter the migrated
   library, so linking the 2stage tests fails with
   `undefined reference to iparam2stage_` etc.

2. **`tests/lapack/reflapack/CMakeLists.txt` glob is also missing
   `*.F`.** Even if the migrator side is fixed, the quad-promoted
   reference library has the same gap (line ~42, `file(GLOB _reflap_f
   ... *.f)` and `*.F90` but no `*.F`).

3. **Multifloats migration of `.F` files produces invalid Fortran.**
   In `/tmp/stg-multifloats/lapack/src/tsytrd_sb2st.F` the rewritten
   PARAMETER assignment (`RZERO = real64x2(limbs=[0.0D0, 0.0_8])`) is
   placed *above* the `IMPLICIT NONE` line. gfortran rejects this with
   "Unexpected IMPLICIT NONE statement". The placement order needs to
   keep the converted-PARAMETER assignments below `IMPLICIT NONE` (and
   below the variable declarations), the way it does for `.f` files.

4. **Reference quad chain segfaults at runtime even when 2stage links.**
   When the manifest+glob fixes were both applied, `dsyev_2stage` and
   `zheev_2stage` segfaulted inside `__GI___libc_free` — heap
   corruption from the `-freal-8-real-16`-promoted Netlib path. The
   D&C variants (`dsyevd_2stage`, `zheevd_2stage`) didn't segfault but
   couldn't be exercised in isolation while sym/Hermitian variants
   crash. Likely an OpenMP-guard or workspace-size mismatch under quad
   promotion.

Every `*_2stage` routine in `doc/PROCEDURES.md` (~12 entries) is
blocked by some combination of the above. Until they're addressed,
skip 2stage routines when picking new test phases.

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
