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
