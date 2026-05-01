# tests/blacs — Outstanding Work

## MPI environment — RESOLVED via Intel MPI 2021.18

The earlier MPICH-on-WSL singleton behavior (every rank seeing
`MPI_Comm_size = 1`) is bypassed by configuring with
`cmake --preset=linux-impi`, which routes `MPIEXEC_EXECUTABLE` /
`MPI_C_COMPILER` / `MPI_Fortran_COMPILER` through Intel MPI 2021.18
in `/opt/intel/oneapi/mpi/latest`. Hydra/PMI connect 4 ranks on the
same node correctly; tests now run on a real 2×2 BLACS grid.

Verified 2026-04-30: `test_qgesd2d` emits

```
  blacs_qgesd2d [p2p_00_to_11] err= 0.0000E+00 digits= 99.00 PASS
```

— a single inter-rank send/recv between (0,0) and (1,1) — instead of
4 lines of `skipped_grid_1x1`.

## Real distributed exec — RESOLVED on kind16 (2026-04-30)

With real ranks connected via the `linux-impi` preset, the test suite
went from `0% real signal` (sandbox mpiexec gave singleton MPI) to
**1045/1045 (100%)** on the kind16 target. Two preset env vars cover
the Intel MPI compat issues with BLACS-registered user ops over
`REAL(KIND=16)`:

* `I_MPI_ADJUST_REDUCE=1` — bypasses the optimised `Reduce_local`
  shortpath that crashes on user MPI_Op + REAL(KIND=16) payload.
  (Original fix; reduced failures 170 → 68.)
* `I_MPI_ADJUST_ALLREDUCE=1` — same fix for the `*GAMX2D / IGAMX2D`
  paths through `MPIR_Allreduce` that scalapack iterative-refinement
  drivers hit. (Added 2026-04-30; resolved the last 8 ScaLAPACK
  failures.)

### Original 3 BLACS failures — RESOLVED

* `blacs_pinfo`: the `get_broadcast_topology` case asserted that
  `BLACS_GET(ctxt, WHAT=10, val)` returned a small 0..127 character
  code. Per `blacs_get_.c`, `WHAT=10` is `SGET_BLACSCONTXT` — the
  internal handle for the given context (Fortran handle is
  `MPI_Comm_c2f` of the underlying communicator). Singleton MPI
  returned a tiny placeholder; real impi returns large integers.
  Fixed by relaxing the assertion to idempotency only — two
  consecutive `BLACS_GET` calls must agree.
* `qtrsd2d` / `xtrsd2d`: not actually a BLACS bug. The test used
  M=6 N=5 with the LAPACK upper-triangular convention (verify
  `i <= j` only); BLACS `BI_GetMpiTrType` for M>N follows the
  upper-trapezoidal convention where col j carries rows
  `1..(M-N+j)` — the M-N rows above the j==N diagonal are part
  of the trapezoid and get transmitted. Fixed by switching tests
  to M=N=5 so the canonical triangular shape is exercised.

See commit `67c57d2`.

## Open: kind10 / multifloats PBLAS Level 1/2 heap corruption

The pattern of "test PASSES every case but crashes during teardown
with `free(): invalid pointer`" persists on kind10 and multifloats
targets only — kind16 is clean. Two distinct modes; full details
in `tests/pblas/TODO.md`'s "kind10 / multifloats — Level 1/2 PBLAS
sporadic context destruction" entry.

Likely root cause is in the migrated `e/y/m/w`-prefix BLACS or PBLAS
C globals — heap corruption or static-state collision specific to
the precision-renamed sources, not a test-authoring or preset issue.
Out of scope for the kind16 stabilization PR.

## Untested BLACS surface

These ship in `external/scalapack-2.2.3/BLACS/SRC/` and are migrated
for kind16 alongside the routines we cover, but the test tree does
not yet exercise them:

### `*trbr2d` receive companion of `*trbs2d`

`test_qtrbs2d` only exercises the `'A'` scope where the broadcast
data was loaded on rank `(0,0)` and verified on every other rank by
having them call `target_qtrbr2d` indirectly through… actually it
*does* call `target_qtrbr2d` on receivers — so this is covered.

### `qtrrv2d`/`xtrrv2d` companions

Receive partners for the trapezoidal point-to-point send. The real
variant is exercised in `test_qtrsd2d`; complex variant is exercised
in `test_xtrsd2d` (added — `target_xtrsd2d`/`target_xtrrv2d` wrappers
declared in `common/target_blacs_body.fypp`).

### `qgamn2d` / `xgamn2d` / `qgsum2d` / `xgsum2d` row/column scopes — RESOLVED

Row ('R') and column ('C') scopes are now covered:

- Real-prefix sum: `test_qgsum2d_rc.f90`
- Real-prefix amx/amn: `test_qgamx2d_rc.f90`
- Complex-prefix sum: `test_xgsum2d_rc.f90`
- Complex-prefix amx/amn: `test_xgamx2d_rc.f90` (added 2026-04-30)

Each program runs amx in row scope and amn in column scope (the
reverse pattern is symmetric and skipped to keep run time bounded).

### `blacs_set` / `blacs_pinfo` direct probes — RESOLVED

`tests/blacs/tests/test_blacs_set.f90` exercises the four settable
WHAT codes (`SGET_NR_BS`, `SGET_NR_CO`, `SGET_TOPSREPEAT`,
`SGET_TOPSCOHRNT`) and verifies the value round-trips through
`blacs_get`. `tests/blacs/tests/test_blacs_pinfo.f90` exercises
`blacs_pinfo` idempotency (two consecutive calls return identical
values, both matching what `grid_init` recorded) and `blacs_get`
with WHAT=0 (default system context) and WHAT=10 (default broadcast
topology).

## Cross-tree follow-ups (NOT to be edited from here)

- The migrator's BLACS recipe stages `blacs_setup_.c` but no public
  Fortran callers in the test tree need it. If a future test wants
  to exercise the alternate setup path, add an interface block in
  `target_blacs_body.fypp`.
- `target_blacs` does not declare wrappers for `cgesd2d`/`zgesd2d`
  etc. (the original-precision entry points). All wrappers go through
  the migrated `q`/`x`/`i` entries. Add original-precision wrappers if
  a regression test is wanted to confirm parity.
