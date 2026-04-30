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

## Auxiliary PBLAS coverage — RESOLVED for abs-matvec family

The 6 auxiliary mat-vec norm entry points (`p[dz]agemv`, `pdasymv`,
`pzahemv`, `p[dz]atrmv`) gain test programs at
`tests/pblas/level2/test_p[dz](agemv|a(sy|he)mv|atrmv).f90` plus
matching wrappers in `common/target_pblas_body.fypp` (interfaces
+ subroutines). References are element-wise hand-coded at quad:
`Y_i := |alpha| * sum_j |op(A)_{i,j}| * |X_j| + |beta * Y_i|`. For
the symmetric / Hermitian / triangular variants the indexing
respects UPLO (with the off-triangle mirrored under abs since
`|conjg(z)| = |z|`); pdatrmv / pzatrmv additionally handle DIAG='U'
by treating on-diagonal entries as 1.

Build verified across all three targets (kind10 / kind16 /
multifloats). Runtime exercise still requires real MPI — same
sandbox limitation that affects the existing pblas tests.

The remaining 7 auxiliary entry points still gap:

- `pdgeadd` / `pzgeadd`: alpha*A + beta*C — reference is element-wise
  add at quad precision.
- `pdtran` / `pztranc` / `pztranu`: transpose (with/without
  conjugation) — reference is a Fortran `transpose()` (real) or the
  same with `conjg` over each element (complex).
- `pdtradd` / `pztradd`: trapezoidal add — reference is element-wise
  add over the upper or lower triangle.

These would follow the same pattern; left as a follow-up.
