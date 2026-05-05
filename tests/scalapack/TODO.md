# tests/scalapack — TODO

Resolved items have moved to `CHANGELOG.md`.

## pzunmrz — SIDE='L' deferred (complex only)

The real-precision SIDE='L' path now passes after four upstream-bug
fixes landed in `recipes/scalapack/source_overrides/`: the
`p[dz]larzb.f` PBxTRAN buffer mismatch, the `p[dz]ormrz.f` post-loop
guard, the `p[dz]larz.f` MPV/NQV undersizing + ZAXPY-stride fix
(carried in from the upstream bug-fix branches; see CHANGELOG and
`doc/UPSTREAM_BUGS.md`), and the analogous `pzlarzc.f` fix. With these
in place, **`test_pdormrz` covers SIDE in {L, R} × TRANS in {N, T}**
(4/4 PASS on kind16 / 2×2 grid).

**`test_pzunmrz` SIDE='L' still fails** by a residual factor of
~1.3-1.8x — distinct from the original ~2.5x. The failure reproduces
on 1, 2, and 4 ranks (so it is not an MPI-distribution issue), and on
both TRANS='N' (1.83x; routes through PZLARZ) and TRANS='C' (1.33x;
routes through PZLARZC). Since the same algorithm structurally works
for the real variant, the remaining bug is complex-specific somewhere
in the PZUNMR3 + PZLARZ / PZLARZC chain — most likely a missed
conjugation or complex-only divergence not yet isolated.

`tests/scalapack/factorization/test_pzunmrz.f90` currently restricts
to SIDE='R' only. Re-add SIDE='L' once the complex-only bug is
identified.
