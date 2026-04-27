# tests/lapack — TODO

Issues observed while adding LAPACK tests that need fixes outside
`tests/lapack/` (migrator, build infra, etc.). Each entry names the
specific routines blocked, the file(s) involved, and what would have
to change.

## Phases L2 / L3 — \*geevx / \*ggev3 / \*ggevx family aborts

The d/z `*geevx`, `*ggev3`, and `*ggevx` tests all abort with
`free(): invalid pointer` shortly after the workspace-correct call.
Pattern is consistent across Phase L2 (L2 surfaced d/zgeevx) and
Phase L3 (L3 added d/zggev3 and d/zggevx, same crash). Workspace
queries return sensible LWORK; the crash happens during the second
call (with allocated WORK) or just after deallocation. Likely a
shared workspace-corruption bug in the migrated `q*ev[x]/q*v3`
chain or in the wrapper itself.

To reproduce: stage kind16 and run any of the failing tests
directly — the abort happens immediately on the first iteration.

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
