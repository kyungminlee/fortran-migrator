# tests/pblas — follow-ups outside this subtree

This file tracks items that were *deferred* during the test-coverage
expansion because they live outside `tests/pblas/`. Each entry names
the file/concern and what should change after the test PR lands.

## doc/PROCEDURES.md — coverage footnote — RESOLVED

`doc/PROCEDURES.md`'s PBLAS section now carries a coverage footnote
listing the 55 of 61 rows exercised by `tests/pblas/`. Counts updated
2026-04-30 from the original 48 (the auxiliary `geadd`, `tradd`,
`tran`, `tranc`, `tranu` cases gained tests since the TODO was
written). Remaining uncovered: 6 auxiliary mat-vec norm stems
(`agemv`, `ahemv`, `asymv`, `atrmv`) — they have no canonical BLAS
analogue and would each need a hand-coded `sum(abs(...))` reference;
see "Auxiliary PBLAS coverage" below.

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
