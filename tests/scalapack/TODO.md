# tests/scalapack — TODO

Resolved items have moved to `CHANGELOG.md`.

## pdormrz / pzunmrz — SIDE='L' deferred

SIDE='R' (both TRANS='N' and TRANS='C/T') passes to target precision on
kind16 / 2×2 grid after two upstream-bug fixes landed in
`recipes/scalapack/source_overrides/p[dz]larzb.f` and
`recipes/scalapack/source_overrides/p[dz]ormrz.f` (see CHANGELOG entry
for the SIDE='R' resolution and `doc/UPSTREAM_BUGS.md` for both bugs).

**SIDE='L' (both TRANS='N' and TRANS='T')** still fails by a factor of
~2.5. The failure reproduces with `K=mA=4` / `mC=nA=K+L=8` (single-block
path, where PDLARZB is *not* called and only the post-loop PDORMR3
fires), ruling out PDLARZB. The bug appears to live in the
PDORMR3 + PDLARZ chain for SIDE='L' (the SIDE='R' chain works
correctly).

Test drivers `tests/scalapack/factorization/test_p[dz]ormrz.f90`
currently restrict to SIDE='R' only. Re-add SIDE='L' once the
underlying PDORMR3 / PDLARZ bug is identified.
