# tests/lapack — TODO

Resolved items have moved to `CHANGELOG.md`. Remaining work: two parked
families with explicit reopen conditions.

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
