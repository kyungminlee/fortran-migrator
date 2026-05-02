# tests/mumps — parallel work groups (2026-05-02)

Source: `tests/mumps/TODO.md` (980 lines, most entries RESOLVED).
This document partitions the surviving open work into five groups that
can be assigned to independent workers. Groups touch disjoint files and
have no implementation-level dependency on each other.

## Open items inventory

After filtering RESOLVED entries from `tests/mumps/TODO.md`:

| ID    | Summary                                                              | Status                          |
| ----- | -------------------------------------------------------------------- | ------------------------------- |
| UNK-3 | Verify `NO_SAVE_RESTORE` setting in staged sources                   | inspection only                 |
| UNK-4 | Verify link-time satisfaction (BLAS/LAPACK/ScaLAPACK/MPI/PORD)       | likely already true             |
| B3    | Full libmpiseq sequential link (qmumps needs Q/X-prefixed P*GET/POT) | partial — alt path open         |
| B4    | `INFOG(20)` cross-target sanity check                                | mmumps build resolved; check open|
| G2    | `target_kind10/`, `target_multifloats/` placeholder dirs             | effectively done — close entry  |
| D1    | Restore error tests via input-sanitising wrapper                     | tests dropped; wrapper missing  |
| B9b   | Test-harness rework for non-kind16 targets                           | three sub-tasks open            |

## Independence map

```
                ┌────────────────────────────────────────┐
                │  tests/mumps/                          │
                │                                        │
   Group A ──▶  │  fortran/   c/   target_multifloats/   │
   (B9b)        │  c/include/test_real_compat.h          │
                │  CMakeLists.txt:52-59 (skip guard)     │
                │                                        │
   Group B ──▶  │  fortran/test_dmumps_infog20.f90 (new) │
   (B4 check)   │                                        │
                │                                        │
   Group D ──▶  │  common/mumps_input_check.f90 (new)    │
   (D1)         │  fortran/test_*mumps_errors.f90 (new)  │
                │                                        │
   Group E ──▶  │  TODO.md (verifications + cleanup)     │
   (UNK-3/4,G2) │                                        │
                └────────────────────────────────────────┘

   Group C ──▶  cmake/CMakeLists.txt section 10
   (B3)         OR external/MUMPS_5.8.2/libseq/mpi.f
                (no test-side changes — mpiexec -n 1 keeps working)
```

No file overlap. Groups A–E proceed in parallel.

---

## Group A — Multifloats test harness (B9b)

**Goal.** Bring multifloats target from 0 mumps tests → full coverage
parity with kind16 / kind10 (currently 23 / 21 tests respectively).

**Files touched.**
- `tests/mumps/c/include/test_real_compat.h` — add `mumps_float64x2`
  POD-struct mode
- `tests/mumps/target_multifloats/c/include/` — already exists
  (`mmumps_c.h`, `wmumps_c.h`, `mumps_c_types.h`)
- `tests/mumps/fortran/test_*mumps_*.f90` — add per-target host↔device
  shim around `id%A` / `id%RHS` assignment + `INFOG/RINFOG/CNTL`
  extraction
- `tests/mumps/c/test_*mumps_c_*.c` — refactor to use POD-aware shim
- `tests/mumps/CMakeLists.txt:52-59` — remove
  `_MUMPS_TESTS_SKIP_EXECUTABLES` guard

**Sub-tasks (loosely coupled, mostly serial within the group).**

### A1 — Multifloats Fortran test executables
Author the host↔device shim — multifloats `id%A` / `id%RHS` are
`TYPE(real64x2)`, not `real(16)`. There is no implicit conversion, so
each assignment site needs `real64x2(real(quad_value, dp), 0.0_dp)` or
similar. Reverse extraction for `INFOG`/`RINFOG`/`CNTL` needs
`real(rinfo_field%hi, 16) + real(rinfo_field%lo, 16)`.

Apply across all 10 driver files in `tests/mumps/fortran/`:
- `test_dmumps_basic.f90`, `test_dmumps_dist_input.f90`,
  `test_dmumps_icntl_io.f90`, `test_dmumps_iref_errchk.f90`,
  `test_dmumps_jobs.f90`, `test_dmumps_jobs_extra.f90`,
  `test_dmumps_nrhs.f90`, `test_dmumps_orderings.f90`,
  `test_dmumps_sym.f90`
- and z-counterparts

Likely cleanest as fypp expansion in `target_multifloats/target_mumps.fypp`
(parallel to `target_kind16/target_mumps.fypp`).

### A2 — Multifloats C-side test harness
Extend `tests/mumps/c/include/test_real_compat.h` with a POD-struct
mode triggered by `TARGET_NAME=multifloats`:
- `TR_LIT(x)` → `((mumps_float64x2){.r=(x), .i=0.0})`
- `TR_FABS(x)` → helper that compares `|x.r| + |x.i|`
- `TR_CMP(a,b)` → componentwise compare
- math: hand-roll the few primitives the tests need (sqrt is unused;
  abs/cmp/literal are sufficient based on
  `test_dmumps_c_basic.c`/`_sym.c` content)

Port the 3 C drivers (`test_dmumps_c_basic.c`, `test_dmumps_c_sym.c`,
`test_zmumps_c_basic.c`) so they compile against the multifloats POD
type.

### A3 — Fortran parity tests portable rewrite
`test_dmumps_c_parity.f90` / `test_zmumps_c_parity.f90` currently embed
`real(16)` in the `bind(C)` interface — kind16-only. Rewrite to use
the per-target precision parameter from `target_mumps.fypp`. For
multifloats this means C-side parity helpers receive a `real64x2`
struct rather than `__float128`; coordinate the bind-C signature with
the C side of `c_parity_helpers.c`.

A3 may share helpers with A1; sketch the shim interface in A1 first.

**Done when:** `ctest -R '^mumps_'` reports the same test count for
multifloats as for kind16.

---

## Group B — `INFOG(20)` cross-target sanity check (B4)

**Goal.** Detect silent mis-sizing of byte-accounting constants in the
`mumps_memory_mod_*.F` overrides — the failure mode B4 warns about that
produces wrong statistics records (not factorisation errors), so it
won't show up in any existing test.

**Files touched.**
- New: `tests/mumps/fortran/test_dmumps_infog20.f90` (or extend
  `test_dmumps_jobs_extra.f90`)

**Implementation.**
1. Factor a fixed problem (n=32, `SYM=0`, `JOB=6`) at the active target
2. Read `INFOG(20)` (peak number of real-precision words used)
3. Divide by `STORAGE_SIZE(real(0,wp))/8` to get a structural element
   count
4. Compare against a reference structural count derived from a kind16
   run baked into the test as a constant
5. Assert agreement within ±5%

A single test program runs at whichever target the build selects;
ctest aggregates across runs.

**Done when:** the test passes for kind10, kind16, and multifloats once
Group A unblocks multifloats execution.

---

## Group C — Full libmpiseq sequential link (B3)

**Goal.** Make qmumps linkable against libmpiseq so `mpiexec` is no
longer required at runtime. Today, the libmpiseq archive itself builds
clean (B3 partial resolution) but linking the migrated qmumps against
it fails — qmumps wants Q-prefixed `pqgetrs_` / `pqgetrf_` / `pqpotrf_`
/ `pqpotrs_` while libmpiseq's `mpi.f` only supplies the D-prefixed
upstream names.

**Files touched.**
- Either `external/MUMPS_5.8.2/libseq/mpi.f` (path C1)
- Or `cmake/CMakeLists.txt` + recipe stripping for libqscalapack /
  libqblacs (path C2)

No test-side changes — current `mpiexec -n 1` path remains the default.
This group runs orthogonally to Group A.

**Two alternative paths — pick one.**

### C1 — Extend `mpi.f` with Q/X analogues (smaller scope)
Add the four routines × precision targets (kind10/kind16/multifloats →
12 Fortran subroutines total) to the libmpiseq stub, each
forwarding/no-op'ing the same way the existing D-prefixed stubs do.

Recommended path — contained to one file, no recipe surgery.

### C2 — Build stripped libqscalapack / libqblacs variants
Add a recipe that produces alternate archives (e.g.
`libqscalapack_seq.a`) omitting the routines libmpiseq already covers:
`descinit_`, `numroc_`, `infog2l_`, `indxg2p_`, `chk1mat_`,
`pchk2mat_`, `pxerbla_`, `descset_`, `blacs_gridinit_`,
`blacs_gridinfo_`, `blacs_gridexit_`. Then link qmumps against the
stripped variant + libmpiseq.

Larger scope; only justify if upstream divergence in C1 becomes
painful.

**Done when:** A test program (one is enough) links and runs without
`mpiexec` and produces the same factorisation as the `mpiexec -n 1`
variant.

---

## Group D — Input-validation wrapper + restore D1 error tests

**Goal.** Re-introduce error-path coverage that was dropped because
MUMPS 5.8.2 doesn't validate negative `N`, out-of-range `IRN`/`JCN`,
or mismatched `NNZ` (it SIGSEGVs instead of returning `INFOG(1) < 0`).

**Files touched.**
- New: `tests/mumps/common/mumps_input_check.f90` — sanitising wrapper
- New: `tests/mumps/fortran/test_dmumps_errors.f90` — restored
- New: `tests/mumps/fortran/test_zmumps_errors.f90` — restored

**Wrapper contract.**
- `id%N >= 1`
- `1 <= IRN(i) <= id%N` for all i in `[1, id%NNZ]`
- `1 <= JCN(i) <= id%N` for all i in `[1, id%NNZ]`
- `size(id%IRN) == id%NNZ`, `size(id%JCN) == id%NNZ`,
  `size(id%A) == id%NNZ`
- Returns a wrapper-defined error code; tests assert against it. Do NOT
  set `id%INFOG(1)` directly — the wrapper layers above MUMPS, doesn't
  pretend to be MUMPS.

**Done when:** `test_*mumps_errors.f90` exercise each invariant and
pass via the wrapper without crashing.

---

## Group E — Verifications + TODO hygiene (UNK-3, UNK-4, G2)

**Goal.** Close out three small bookkeeping items in `TODO.md`.

**Files touched.**
- `tests/mumps/TODO.md` only

**Tasks.**

### E1 — UNK-3: NO_SAVE_RESTORE inspection
After `pyengine stage --target kind16 --libraries mumps`, grep
`/tmp/stg-q/mumps/src/` for `#if ! defined(NO_SAVE_RESTORE)` and
report which branch is active in the staged tree. Update UNK-3 with
the resolved status and mark it RESOLVED.

### E2 — UNK-4: link-time satisfaction
Confirm the existing build evidence is sufficient: B5 already shows
`libqmumps-gfortran-13.a` at ~22 MB built clean. Spot-check that no
PORD symbols leak as unresolved by running `nm libqmumps-gfortran-13.a
| grep -i pord | grep ' U '` and noting the result. Update UNK-4
with the verified status.

### E3 — G2: close out
Both `target_kind10/target_mumps.fypp` and
`target_multifloats/target_mumps.fypp` already exist. Mark G2
RESOLVED with a one-liner pointing at the on-disk evidence.

**Done when:** all three entries carry a RESOLVED tag with date.
Quickest group — ~30 min sitting.

---

## Estimated workload split

Rough estimates based on file count, novelty, and known unknowns:

| Group | %    | Rationale                                                                  |
| ----- | ---- | -------------------------------------------------------------------------- |
| A     | ~60% | Shim design + 10 Fortran drivers + 3 C drivers + parity rewrite. Largest because every assignment site to `id%A`/`id%RHS`/INFOG needs the host↔device wrap, and the POD-struct C path has no precedent in `test_real_compat.h`. |
| B     | ~7%  | One new test, well-bounded. Hardest part is picking the ±5% threshold; the code itself is short. |
| C     | ~20% | C1 path: 12 forwarding stubs (4 routines × 3 targets) + verify link. Hidden risk of second-order missing symbols once the first batch lands. |
| D     | ~12% | Wrapper (~half day) + two restored error tests (one each). Straightforward but multiplies across two precisions. |
| E     | ~1%  | Pure inspection + TODO edits. ~30 min sitting.                              |

**Variance.** Group A could push to 65–70% if the multifloats Fortran
shim turns out to need per-field handling for `RINFO`/`CNTL` arrays
(not just scalars). Group C could push toward 25% if C1 surfaces
duplicate-symbol issues — the open-part of B3 is flagged "non-trivial"
in the TODO.

## Suggested assignment

| Worker | Group | Notes                                                    |
| ------ | ----- | -------------------------------------------------------- |
| 1      | A     | Largest; A1/A2/A3 mostly serial within the group         |
| 2      | B     | Small additive test                                      |
| 3      | C     | Build-system; recommend C1 (smaller scope)               |
| 4      | D     | Wrapper + 2 restored error-path tests                    |
| 5      | E     | Verifications + TODO cleanup; quickest (~30 min)         |

Coordination is minimal: Group A's removal of the
`_MUMPS_TESTS_SKIP_EXECUTABLES` guard in `CMakeLists.txt:52-59` is the
only change another group could conflict with. Group B should target
kind10/kind16 first and pick up multifloats once Group A lands.
