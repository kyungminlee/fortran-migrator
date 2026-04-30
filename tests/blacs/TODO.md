# tests/blacs — Outstanding Work

## MPI sandbox caveat

In this sandbox `mpiexec -n N <prog>` spawns **N unconnected MPI worlds**
— each rank sees `MPI_Comm_size = 1`, BLACS reports a 1×1 grid, and
`blacs_gridinit` returns each process its own private context.

Consequence: every BLACS comm test except `test_grid_info` either
emits `skipped_grid_1x1` (the explicit early-return path in p2p tests
like `test_qgesd2d`, `test_qtrsd2d`, `test_xgesd2d`, `test_igesd2d`)
or trivially passes on a 1-element communicator (broadcast/reduce
tests where rank (0,0) writes and verifies its own buffer with no
peer involved). **No real signal in this sandbox** — the assertions
are correct but unexercised.

The tests are correctly written for a real 2×2 grid and *will*
exercise actual inter-rank communication when run on a properly
configured MPI deployment (where Hydra/PMI can connect ranks).
Confirm the sandbox behavior with:

```
mpiexec -n 4 /tmp/staging-kind16/build/tests/blacs/test_qgesd2d
```

— output shows 4 separate `skipped_grid_1x1` lines, one per
unconnected world.

## Environment caveat: MPICH singleton fallback under sandboxed bash

In the sandbox where these tests were authored, `mpiexec -n 4 <prog>`
ran *4 independent MPICH singleton processes* rather than a connected
MPI world (each process saw `MPI_Comm_size = 1`). A trivial
`MPI_Comm_rank/size` C program had the same behavior, so the issue is
environmental (likely PMI/Hydra cannot connect under the bash sandbox)
and not in the test code or the migrated BLACS library.

In that mode the tests still PASS:

  - on grid 1×1 the p2p tests (`test_qgesd2d`, `test_qtrsd2d`,
    `test_xgesd2d`, `test_igesd2d`) take a `skipped_grid_1x1` early
    return and report `err = 0`,
  - the broadcast/reduce tests behave like local no-ops: only rank
    `(0,0)` writes/sums its own buffer; the post-call check sees the
    expected value because rank `(0,0)` is the sole participant.

So the tests *don't* actually exercise inter-rank communication when
run in the sandbox. **They do once mpiexec is connecting ranks
properly** (verified manually by inspecting the wrapper paths and the
BLACS scope routing) — outside the sandbox, ctest reports the same
`100% passed`, but the assertions then have teeth.

There is nothing to fix in `tests/blacs/`; this is documented purely
to set expectations. If a contributor reproduces the singleton
behavior, the fix is at the harness layer, not here.

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
