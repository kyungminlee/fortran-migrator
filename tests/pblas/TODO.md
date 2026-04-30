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

The remaining 7 auxiliary entry points are also covered (TODO entry
was stale through a prior edit; verified 2026-04-30):

- `pdgeadd` / `pzgeadd`: `tests/pblas/level3/test_p[dz]geadd.f90`
  (element-wise `beta*C + alpha*op(A)` reference at quad).
- `pdtran` / `pztranc` / `pztranu`: `tests/pblas/level3/test_p[dz]tran[cu]?.f90`
  (Fortran `transpose()` / `transpose(conjg(...))` references at quad).
- `pdtradd` / `pztradd`: `tests/pblas/level3/test_p[dz]tradd.f90`
  (UPLO-restricted element-wise add at quad).

Wrappers in `common/target_pblas_body.fypp` (`target_pdgeadd`,
`target_pdtran`, `target_pdtradd`, `target_pzgeadd`, `target_pztranc`,
`target_pztranu`, `target_pztradd`) and a corresponding test file each.
Build verified across kind10 / kind16 / multifloats; runtime exercise
under real distributed MPI lands once `cmake/CMakePresets.json`'s
`linux-impi` preset is invoked.

## Real-MPI exercise via `linux-impi` preset — 2026-04-30

`cmake --preset=linux-impi` routes Intel MPI 2021.18 from
`/opt/intel/oneapi/mpi/latest`; `MPIEXEC_PREFLAGS` injects
`-genv I_MPI_ADJUST_REDUCE=1` to skip the impi shortpath that
mishandles BLACS-registered user ops over REAL(KIND=16). With that:

- Most pblas tests PASS every documented case under a real 2×2 grid.
- 41 still hit `free(): invalid pointer` during teardown OR crash
  mid-suite on specific transpose paths (e.g. `pdgemm` case 8 with
  TRANSA=N, TRANSB=T).

The teardown crashes are BLACS/PBLAS-level MPI lifecycle bugs —
likely `MPI_Op_free` / `MPI_Type_free` ordering vs `BLACS_EXIT(0)` /
`MPI_Finalize`. Tracked in `tests/blacs/TODO.md`; not test bugs.
