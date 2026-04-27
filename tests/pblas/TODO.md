# tests/pblas — follow-ups outside this subtree

This file tracks items that were *deferred* during the test-coverage
expansion because they live outside `tests/pblas/`. Each entry names
the file/concern and what should change after the test PR lands.

## doc/PROCEDURES.md — coverage column

`doc/PROCEDURES.md` enumerates all 61 PBLAS (family, stem) rows with
their migrated names per target, but has no column indicating which
ones are exercised by `tests/pblas/`. With this PR, 48 of the 61
rows are exercised on the kind10 and kind16 targets:

- All 7 Real-prefix Level 1 stems
- All 4 Complex-prefix Level 1 stems
- Both Mixed real-from-complex stems (asum, nrm2)
- The Complex-vector with real-scalar stem (zdscal)
- All 7 Real-prefix Level 2 stems (gemv, ger, symv, syr, syr2, trmv, trsv)
- All 8 Complex-prefix Level 2 stems (gemv, hemv, gerc, geru, her, her2, trmv, trsv)
- All 6 Real-prefix Level 3 stems (gemm, symm, syrk, syr2k, trmm, trsm)
- All 8 Complex-prefix Level 3 stems (gemm, hemm, symm, herk, her2k, syrk, syr2k, trmm, trsm)

Rows *not* exercised — auxiliary PBLAS routines without a direct
serial-BLAS analogue:

- `agemv`, `ahemv`, `asymv`, `atrmv` (auxiliary mat-vec norms)
- `geadd`, `tradd` (general / triangular add)
- `tran`, `tranc`, `tranu` (transpose / conjugate-transpose)

A coverage column or an explicit "tested-by" footnote would document
this without forcing readers to grep `tests/pblas/level*/`.

## Multifloats `extern "C"` linkage

`tests/pblas/README.md` (Open Issues) and `recipes/pblas.yaml`
(`header_patches: when: multifloats`) already document the situation:
the migrated PBLAS .c sources compile as C++ for the multifloats
target, so the Fortran-callable entry points are name-mangled and
unreachable from `target_pblas.f90`. With 28 new tests, the gap is
larger (40+ tests skipped instead of 20) and the cost-benefit of
fixing the linkage has shifted. Suggested follow-up:

- Confirm the `extern "C"` wrapping in `recipes/pblas.yaml`'s
  `header_patches` and `c_migrator_overrides` actually round-trips
  through to the per-file `<name>_.c` definitions.
- If not, extend the C migrator's body wrap (or add a `body_wrap`
  field to recipes) so that `extern "C" { … }` brackets every entry
  point alongside its declaration.

Once that is in place, `tests/pblas/CMakeLists.txt` can drop the
"return early when `C_AS_CXX` is on" guard and multifloats joins
kind10/kind16 in the precision report.

## C migrator: cmplx16 alias not applied to local-variable declarations

`recipes/pblas.yaml` declares

```yaml
c_type_aliases:
  - from: [cmplx, cmplx16]
    to: 'cmplx{RPU}'
```

so that `cmplx16` (Netlib's `typedef double cmplx16[2]`) gets renamed
to the kind-specific complex typedef (`cmplxQ` for kind16,
`cmplxE` for kind10, `cmplxDD` for multifloats). The alias is applied
to function-signature parameter types, but **not** to local-variable
declarations inside function bodies. Consequence: the migrated
`PB_Ctzhemm.c` and `PB_Ctzher2k.c` still declare

```c
cmplx16        Calph16;
```

which is a 16-byte buffer (2 × `double`) on the kind16 target, even
though the surrounding code passes `Calph16` to the migrated `gemm`
expecting a 32-byte `cmplxQ`. The buffer is also written to by
`PB_Cconjg`, whose body uses hardcoded `(double*)` casts and was not
upgraded by the migrator either:

```c
case DCPLX:
    ((double*)(CALPHA))[REAL_PART] =  ((double*)(ALPHA))[REAL_PART];
    ((double*)(CALPHA))[IMAG_PART] = -((double*)(ALPHA))[IMAG_PART];
```

For the kind16 target the type marker is still `DCPLX` (set in
`pb_cxtypeset.c`), so this case is taken — but it copies only 16
bytes from `ALPHA`, treating the first 16-byte half of a `cmplxQ` as
two adjacent `double`s. The result `Calph16` contains a corrupted
half-precision reconstruction of `conjg(alpha)`, and any kernel call
that uses `Calph16` (rank-2k, Hermitian matmul, etc.) gets garbage
for `conjg(alpha)`.

`tests/pblas/level3/test_pzher2k.f90` exercises exactly this path
(complex `alpha` with non-zero imaginary part) and surfaces the bug
as max-rel-err ≈ 1.0+ with tolerance ≈ 1e-30. The test is
intentionally left in place so the failure is visible; once the
migrator bug is fixed, the test should go green automatically.

The test now also exercises the additional UPLO/TRANS combinations
(UPLO ∈ {'U','L'} × TRANS ∈ {'N','C'}, four combos × three sizes =
12 cases). Every combo routes through the same PB_Ctzher2k stack
allocation and PB_Cconjg call site, so all 12 cases also fail with
max-rel-err ≈ 1.0+ on kind16. They are kept so the migrator fix can
be validated across PB_Cptzher2k{LU,LL,UN,UC,LN,LC} simultaneously.

Fix should land in two places, both outside `tests/pblas/`:

1. `src/pyengine/c_migrator.py` — extend `c_type_aliases` to rewrite
   local-variable declarations, not just function-signature types.
2. The same migrator pass should also rewrite the hardcoded
   `(double*)` / `(float*)` casts inside `PB_Cconjg.c` (and any
   sibling helpers) to the kind-correct real-pointer type.

Once the fix lands, also re-check whether `PB_Cpsymm{AB,BC}.c`,
`PB_Cptrmm{AB,BC}.c` and friends — which similarly stack-allocate
`cmplx16 talpha[…]` — exercise their conjg paths during existing
tests, and add a small smoke test if the path is otherwise unreached.

## Follow-up: re-run pzher2k after migrator fix

Once the `cmplx16` aliasing / `PB_Cconjg` cast fix above lands in
`src/pyengine/c_migrator.py`, re-run

```bash
ctest --test-dir <build> -R '^pblas_test_pzher2k$' --output-on-failure
```

against a fresh kind16 staging build. Expected outcome: max-rel-err
drops from ≈ 1 to ≤ `64 * 8 * k * eps` (≈ 1e-30 for the largest case)
and the test goes green. If it still fails, the fix is incomplete —
check the kind10 and multifloats targets too, since the alias rule
covers all three, and remove the explanatory header comment at the
top of `tests/pblas/level3/test_pzher2k.f90` once the bug is gone.

## Auxiliary PBLAS coverage

The 13 auxiliary entry points listed above (`agemv`, `geadd`, `tran`,
…) have no canonical BLAS reference, so they are not part of this
PR's scope. They are still part of the migrated PBLAS surface and
worth covering with bespoke tests:

- `pdgeadd` / `pzgeadd`: alpha*A + beta*C — reference is element-wise
  add at quad precision.
- `pdtran` / `pztranc` / `pztranu`: transpose (with/without
  conjugation) — reference is a Fortran `transpose()` (real) or the
  same with `conjg` over each element (complex).
- `pdtradd` / `pztradd`: trapezoidal add — reference is element-wise
  add over the upper or lower triangle.
- `pdagemv` / `pdasymv` / `pdatrmv` / `pzagemv` / `pzahemv` /
  `pzatrmv`: produce a vector of row/column 1-norms of the
  underlying mat-vec — reference is `sum(abs(...))` over each
  row/column at quad precision.

These tests would fit directly into `tests/pblas/level{2,3}/`
without further infrastructure, but each needs a small inline
reference computation in lieu of a Netlib BLAS analogue. Worth
adding in a follow-up PR.
