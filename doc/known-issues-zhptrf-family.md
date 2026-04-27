# Migrated Hermitian-packed family (`zhptrf` / `zhptrs` / `zhpsv` / `zhptri`) diverges from Netlib reference

**TL;DR.** The migrated `xhptrf` (kind16) and `yhptrf` (kind10) entry
points produce factors / solutions that diverge from the
`-freal-8-real-16` quad-promoted Netlib `zhptrf` by ~1e-7 — orders of
magnitude beyond floating-point roundoff. The same test passes for
`zsptrf` (complex symmetric packed), `zpptrf` (Hermitian PD packed), and
all dense `zhetrf` variants, so the bug is specific to the
Bunch-Kaufman packed Hermitian path. Multifloats happens to pass — its
~32-digit precision masks whatever rounding-sensitive divergence is
present.

## Reproducer

Tests live under `tests/lapack/{factorization,linear_solve}/test_zhptr*.f90`
in branch where they were briefly added; they were dropped from the
test suite to keep CI green pending investigation. The wrappers and
ref interfaces remain in `target_lapack_body.fypp` and
`ref_quad_lapack.f90` so the tests can be reinstated trivially once
the migrated code is fixed.

```bash
# Re-add the four test programs and run:
ctest --test-dir <build> -R 'lapack_test_(zhptrf|zhptrs|zhpsv|zhptri)$' \
      --output-on-failure
```

Observed: max-rel-err ~1e-7 on n=16/32/64 for kind10 and kind16 targets;
multifloats passes at full precision.

## Suspected cause

Likely a mishandled conjugation in the migrated Hermitian packed
routines — the algorithm reads `AP(diag)` as real but writes back via
the standard complex conversion path, which would lose the
"diagonal must be real" invariant in the migrated kind10/kind16 builds
but not in multifloats (where the extra precision absorbs the error).

Needs inspection of the migrated `xhptf2` / `yhptf2` (the unblocked
panel) and any explicit `real()` calls that may have been dropped.
