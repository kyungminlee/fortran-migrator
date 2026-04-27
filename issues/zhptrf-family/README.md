# Migrated Hermitian-packed family (`zhptrf` / `zhptrs` / `zhpsv` / `zhptri`) diverges from Netlib reference

| | |
|---|---|
| **First observed at commit** | `1e13e0b` (branch `tests`, 2026-04-26) |
| **Originally documented at** | `ae39cbf` — Phase 10 commit that dropped these tests |
| **Status when this directory was created** | open / unfixed |

> When the upstream code in `external/lapack-3.12.1/SRC/{zhptrf,zhptrs,zhpsv,zhptri}.f`
> or the migrator's recipe for the Hermitian-packed family changes, re-run the
> staged test programs in `tests-dropped/` against `HEAD`. If they now pass on
> all three targets, this issue is fixed and the directory can be removed.

## TL;DR

The migrated `xhptrf` (kind16) and `yhptrf` (kind10) entry points produce
factors / solutions that diverge from the `-freal-8-real-16` quad-promoted
Netlib `zhptrf` by **~1×10⁻⁷ in absolute relative error**, regardless of `n`.
That is many orders of magnitude beyond floating-point roundoff — it is an
algorithmic divergence, not precision noise. The same test passes for
`zsptrf` (complex symmetric packed), `zpptrf` (Hermitian PD packed), and all
dense `zhetrf` variants, so the bug is specific to the Bunch-Kaufman packed
Hermitian path. Multifloats happens to pass — its ~32-digit precision masks
whatever rounding-sensitive divergence is present.

## Reproducer

```bash
# Stage and build any one target
PYTHONPATH=src python -m pyengine stage /tmp/stg-kind10 \
    --target kind10 --libraries blas blacs lapack
cmake -S /tmp/stg-kind10 -B /tmp/stg-kind10/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg-kind10/build -j$(nproc)

# Drop the four test programs into the staged tree
cp issues/zhptrf-family/tests-dropped/test_zhpt*.f90 \
   /tmp/stg-kind10/tests/lapack/factorization/   # zhptrf, zhptri
cp issues/zhptrf-family/tests-dropped/test_zhp{sv,trs}.f90 \
   /tmp/stg-kind10/tests/lapack/linear_solve/
cmake --build /tmp/stg-kind10/build -j$(nproc)

ctest --test-dir /tmp/stg-kind10/build \
      -R 'lapack_test_(zhptrf|zhptrs|zhpsv|zhptri)$' --output-on-failure
```

Observed at the commit referenced above:

| target | n=16 | n=32 | n=64 |
|---|---|---|---|
| kind10 (`yhptrf`) | 1.99e-7 | 3.39e-7 | 4.17e-7 |
| kind16 (`xhptrf`) | similar magnitude | similar | similar |
| multifloats (`vhptrf`) | passes (≪ tol) | passes | passes |

Tolerance applied: `tol = 16 * n² * target_eps`.

## Why we suspect a real-vs-complex conversion error on the diagonal

- **Specific to Hermitian packed.** `zsptrf` (complex symmetric packed) and
  `zpptrf` (Hermitian PD packed) both pass cleanly; only the Bunch-Kaufman
  Hermitian-packed routines diverge.
- **Multifloats passes; kind16 fails.** kind16 should be bit-identical to the
  reference because both are quad-precision promoted Netlib. The fact that
  it isn't says the migrated code path differs in code, not in arithmetic.
- **Hermitian-packed BK enforces real diagonal.** If the unblocked panel
  inside `xhptrf`/`yhptrf` lost an explicit `real()` projection during
  migration but `zhetf2` (dense, separately tested) kept it, the diagonal
  could drift imaginary, change pivot decisions, and propagate. Multifloats'
  extra precision would dilute the imaginary contamination below the
  pivot-decision threshold and mask the bug.

## Files in this directory

- `netlib/` — pristine Netlib LAPACK 3.12.1 sources used as the reference
  (compiled with `-freal-8-real-16` for the differential test):
  `zhptrf.f`, `zhptrs.f`, `zhpsv.f`, `zhptri.f`
- `kind10/` — migrated 80-bit extended sources (`yhptrf.f`, `yhptrs.f`,
  `yhpsv.f`, `yhptri.f`)
- `kind16/` — migrated `REAL(KIND=16)` sources (`xhptrf.f`, `xhptrs.f`,
  `xhpsv.f`, `xhptri.f`)
- `multifloats/` — migrated `real64x2` sources (`vhptrf.f`, `vhptrs.f`,
  `vhpsv.f`, `vhptri.f`) — these pass and are included as the working
  reference for diff'ing
- `tests-dropped/` — the four test programs that were written, observed to
  fail, and excluded from the test suite at `ae39cbf`. They depend on the
  `target_zhptrf`/`target_zhptrs`/`target_zhpsv`/`target_zhptri` wrappers
  that were kept in `tests/lapack/common/target_lapack_body.fypp`, so they
  drop in cleanly.

## Suggested first diff

```bash
diff issues/zhptrf-family/multifloats/vhptrf.f \
     issues/zhptrf-family/kind16/xhptrf.f
```

If the migration is mechanical, the two files should differ only in
type-name substitutions (`real64x2` ↔ `real(kind=16)`,
`cmplx64x2` ↔ `complex(kind=16)`) and the prefix letter. Any structural
divergence — especially around handling of `AP(diag)` or the BK pivot
selection block — is the place to look.
