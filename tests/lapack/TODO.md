# tests/lapack — TODO

Resolved items have moved to `CHANGELOG.md`. Remaining work: two parked
families with explicit reopen conditions.

## Older blocked-T routines that need a wider scope to test — PARKED

The following routines were tested in Phase 35 / Phase 36, but only the
|R| (or |L|) factor is compared — the lower-triangle reflector storage
and the workspace `T` (whether the reserved workspace of `dgeqr` or the
explicit block-T of `dgeqrt`) differ between implementations because
the recursion / blocking choices in the migrated and reference paths
diverge:

- `dgeqr` / `zgeqr`, `dgelq` / `zgelq` (Phase 35)
- `dgeqrt` / `zgeqrt`, `dgelqt` / `zgelqt` (Phase 36)

For the matching `*gemqr` / `*gemlq` and `*gemqrt` / `*gemlqt` apply
routines the C output is canonical and matches cleanly — those tests
do compare the full result. So the factorization+apply pipeline is
covered end-to-end; only the intermediate reflector/T representations
are loose.

**Reopen condition**: a migrator change (or recipe override) that pins
both paths to the same blocking heuristic. Until then the |R|/|L|
comparison is the strongest assertion the test framework can make.

## sb2st kernels (2 routines) — PARKED

`dsb2st_kernels`, `zhb2st_kernels` — SBR (Successive Band Reduction)
inner kernels, normally driven by `*sb2st` / `*hb2st`. Testing them in
isolation requires constructing a partial band-reduction state that
matches the kernel's mid-stream invariants — non-trivial.

These are the only two user-facing routines without a test driver per
the audit. The orchestrator paths (`*sb2st` / `*hb2st`) ARE tested
end-to-end, and they exercise the kernels on every call, so coverage
is effectively transitive.

**Reopen condition**: a regression or precision divergence
attributable to the kernels that the orchestrator tests fail to catch.
Until then, standalone drivers aren't worth the state-construction
work.
