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
mishandles BLACS-registered user ops over REAL(KIND=16).

Three classes of failures surfaced and were fixed (commit 303df23):

1. **Buffer-size mismatch on non-owning ranks** — `gen_distrib_vector`
   zeroed `loc_n` on `my_col != 0` and allocated `max(1, 0) = 1`,
   but the descriptor's LLD was computed without that gate. Wrappers
   reading `x(1:lld*N)` overran. Fix: allocate full `max(1, full_loc_n)`
   on every rank.

2. **`desc(9)*desc(4)` overshoots actual buffer** — every wrapper
   computed `LLD * full-N`, not `LLD * locn`. Fixed by introducing
   `desc_buf_size(desc)` that calls `BLACS_GRIDINFO + numroc_local`.
   Replaced 76 call sites in `target_pblas_body.fypp`. The naive form
   happened to be correct under singleton MPI (npcol=1 → locn=N).

3. **Cabs1 vs Euclidean magnitude** — references for the auxiliary
   `pz{agemv,ahemv,atrmv}` used Fortran's intrinsic `abs(complex)`
   (returns `sqrt(Re^2+Im^2)`), but PBLAS uses Cabs1 = `|Re|+|Im|`.
   References now compute Cabs1 explicitly. zahemv additionally reads
   only `|Re(A_jj)|` for Hermitian diagonal entries.

## kind10 / multifloats — Level 1/2 PBLAS sporadic context destruction — RESOLVED 2026-05-01

This entry described 18 PBLAS L1/L2 failures on kind10 and 26 on
multifloats (heap corruption at exit, sporadic mid-loop
`Cblacs_gridinfo`-returns-`-1`, etc.). When re-tested 2026-05-01 after
the audit-branch merge (commit `cc54c42`) and the multifloats
migrator polish (`0caddd1`, `bf33cde`, `1a609a0`, `0cfb366`), all
three targets are clean:

* **kind16/impi**: 1045/1045 (100%, incl. mumps tests).
* **kind10/impi**: 1022/1022 (100%, mumps tests excluded — pre-existing
  hardcoded q/x prefixes in `tests/mumps/CMakeLists.txt`).
* **multifloats/impi**: 1022/1022 (100%, mumps excluded — same).

The previously-described L1/L2 failure modes (SIGABRT-on-exit and
sporadic `BI_MyContxts[ctxt]` NULL-flip) no longer reproduce in the
current build. Most likely the audit-branch merge brought changes that
fixed the underlying state collision, but the original symptoms also
manifested as Heisenbugs (timing-dependent across rebuilds), so this
entry is closed conservatively — re-open if the failures recur.

The 3 multifloats `pzdtsv` / `pzdttrf` / `pzdttrs` failures noted
alongside this entry were a separate bug (CMPLX-of-complex identity
mishandling in the migrator) — fixed in `0caddd1` and verified at
~31-digit accuracy.
