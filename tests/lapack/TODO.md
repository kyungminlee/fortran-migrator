# tests/lapack — TODO

Issues observed while adding LAPACK tests that need fixes outside
`tests/lapack/` (migrator, build infra, etc.). Each entry names the
specific routines blocked, the file(s) involved, and what would have
to change.

## Phase L2 — three test-side failures left in place

Three new tests landed RED with Phase L2:

- `test_dgeevx` and `test_zgeevx`: abort with `free(): invalid pointer`
  during the first `dgeevx`/`zgeevx` call. Workspace query path looks
  correct (LWORK from `WOPT(1)`, IWORK sized 2*N-2, LDVL=LDVR=1 when
  JOBVL/JOBVR='N' per LAPACK manual). May surface when
  `BALANC='N' + SENSE='N'` is combined with the migrated path on
  kind16; investigate with valgrind / address sanitizer.

- `test_ztrsen`: the `n=24` case reports `max_rel_err ≈ 0.85` while
  `n=12` and `n=32` pass with `err = 0`. Almost certainly a
  comparison-side bug — sorting moduli can swap selected vs
  unselected eigenvalues with equal magnitudes — but the asymmetry
  vs `n=12/32` is unusual. Recheck with explicit per-position
  comparison after a stable secondary sort.

Tests are left in the build per the failing-test-stays-visible
convention. Fix in a follow-up phase.

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
