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

## Open: Intel MPI compat issues surfaced under real distributed exec

With real ranks connected via the `linux-impi` preset, the test suite
went from `0% real signal` (sandbox mpiexec gave singleton MPI) to
`932/1045 passed` (89%). 113 failures remain, distributed across
BLACS / PBLAS / ScaLAPACK.

`I_MPI_ADJUST_REDUCE=1` (force 'binomial' algorithm; skip the Intel
MPI shortpath that dispatches a user-registered `MPI_Op` through
intrinsic `MPIR_SUM`) is wired into the preset's `MPIEXEC_PREFLAGS` —
that alone reduced failures from 170 to 68.

### Remaining BLACS failures (3 of ~21)

```
blacs_test_blacs_pinfo
blacs_test_qtrsd2d      blacs_test_xtrsd2d
```

- `blacs_pinfo`: probable test-program teardown / two-phase-init
  ordering bug under real impi.
- `qtrsd2d` / `xtrsd2d`: numerical garbage on the trapezoidal p2p
  (`err=************, digits=*****`) — receiver decodes a different
  byte layout than sender wrote. Suspect mismatch in
  `BI_qtrsd2d` / `BI_qtrrv2d` MPI_Datatype packing under impi's send
  path. Needs inspection of the BI_q*.c packers vs Intel MPI's
  contiguous-only shortcut for KIND=16.

### Remaining PBLAS / ScaLAPACK failure pattern

A common pattern: the test program PASSES every documented case but
crashes during teardown (`free(): invalid pointer`) or mid-suite on a
specific transpose path. The correctness signal is correct; the crash
is during MPI/BLACS resource cleanup or between-cases reset.

Likely root cause: `MPI_Op_free` / `MPI_Type_free` ordering, or
`BLACS_EXIT(0)` followed by `MPI_Finalize` double-frees a
communicator that one of the BLACS-registered ops also references.
Intel MPI surfaces this where MPICH's permissive cleanup did not.

This is BLACS/PBLAS-level MPI lifecycle work — not test authoring,
not preset wiring. Tracked here for visibility.

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
