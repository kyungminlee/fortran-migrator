# tests/lapack — TODO

Issues observed while adding LAPACK tests that need fixes outside
`tests/lapack/` (migrator, build infra, etc.). Each entry names the
specific routines blocked, the file(s) involved, and what would have
to change.

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
