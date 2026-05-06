# tests/lapack — TODO

Resolved items have moved to `CHANGELOG.md`. Remaining work: parked
items with explicit reopen conditions.

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

## kind4 / kind8 baseline missing the orbdb3 override — PARKED

Upstream LAPACK 3.12.1 ships `{s,d,c,z}orbdb3.f` with a typo on the
`{S,D,Z,C}ROT` call inside the `I .GT. 1` branch — passes `LDX11` as
the `INCY` for `X21`, where it should be `LDX21`. When `LDX11 /=
LDX21` (the orbdb3 regime, `M-P < P`), `*ROT` strides past `X21` and
writes into the next heap chunk. Manifestation varies with heap
layout: `zunbdb3` SIGSEGVs in the wrapper's auto-deallocate;
`dorbdb3` silently corrupts adjacent `theta` / `phi` chunks and
presents as a "wrong digits" precision result. Confirmed both fail
modes resolve when the override is overlaid into the kind8 baseline's
`_reflapack_src/`.

The migrator already carries the one-character fix in
`recipes/lapack/source_overrides/{s,d,c,z}orbdb3.f` and applies it for
every migrated target (kind10, kind16, multifloats). The kind4 / kind8
**baseline** path (`migrator stage --target kind{4,8}`) skips
migration entirely and copies upstream `external/lapack-3.12.1/SRC/`
verbatim into `_reflapack_src/`, so the std `lapack` archive that the
baseline ctest cycle links against still has the upstream typo.

This is a migrator/staging gap, not a test-side bug — the wrappers
and refblas/reflapack pipelines are correct. Routing here per the
"stay inside tests/lapack/" rule.

**Reopen condition**: when the kind4 / kind8 baseline column is
promoted to part of the published precision matrix. At that point
the staging path needs to overlay `recipes/lapack/source_overrides/`
on top of `_reflapack_src/` after the upstream copytree (and
analogous treatment for any later `recipes/<lib>/source_overrides/`
directories). Until then, the 6-failure baseline noted in
`tests/KIND48_BASELINE_STATUS.md` (2 of which trace to this bug) is
acceptable.
