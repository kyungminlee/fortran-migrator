# tests/mumps — DEFERRED TEST IMPLEMENTATION PLAN

> **Status: NEEDS TO BE REVIEWED AGAIN before implementation begins.**
>
> This file captures the test design from the planning session of 2026-04-29.
> Test sources, helpers, and `CMakeLists.txt` are deferred until the
> build-pipeline integration (pyengine `LIBRARY_ORDER` + `cmake/CMakeLists.txt`)
> is verified end-to-end and the open questions below are resolved.

## Open pre-implementation questions

### UNK-1 — Migrated symbol names — RESOLVED (verified against staged output)

Migrator's `prefix: direct` rule **replaces** the leading precision letter
(D/Z for double-real / double-complex) with the target letter (Q/X for
kind16) for **filenames, modules, and subroutines** — but **leaves
type names unchanged**. Verified by inspecting `/tmp/stg-q/mumps/src/`
after `pyengine stage --target kind16 --libraries mumps`.

| Upstream                  | Migrated kind16        |
| ------------------------- | ---------------------- |
| `DMUMPS` subroutine       | `QMUMPS`               |
| `ZMUMPS` subroutine       | `XMUMPS`               |
| `DMUMPS_STRUC_DEF` module | `QMUMPS_STRUC_DEF`     |
| `ZMUMPS_STRUC_DEF` module | `XMUMPS_STRUC_DEF`     |
| `dmumps_struc.h` header   | `qmumps_struc.h`       |
| `zmumps_struc.h` header   | `xmumps_struc.h`       |
| `DMUMPS_STRUC` type       | **`DMUMPS_STRUC`** (unchanged — see note) |
| `ZMUMPS_STRUC` type       | **`ZMUMPS_STRUC`** (unchanged — see note) |
| `A`, `RHS`, `COLSCA`, `RINFO`, etc. fields | Promoted to `REAL(KIND=16)` / `COMPLEX(KIND=16)` |

**Why the type name doesn't move:** the migrator renames the wrapping
module (`MODULE QMUMPS_STRUC_DEF { INCLUDE 'qmumps_struc.h' }`) and the
header file, but the `TYPE …` inside `qmumps_struc.h` still reads `TYPE
DMUMPS_STRUC`. The wrapping module is enough to disambiguate at use
sites; cf. `qmumps_driver.F:197` — `TYPE (DMUMPS_STRUC), TARGET :: id`.

Test code therefore writes:

```fortran
USE QMUMPS_STRUC_DEF                    ! migrated module
TYPE(DMUMPS_STRUC), TARGET :: id        ! type name unchanged
CALL QMUMPS(id)                         ! migrated subroutine
```

For the complex side: `USE XMUMPS_STRUC_DEF`, `TYPE(ZMUMPS_STRUC) :: id`,
`CALL XMUMPS(id)`.

Same pattern applies to other arithmetic-lettered MUMPS symbols the tests
will touch (e.g. `QMUMPS_INTR_TYPES` module wrapping unchanged
`DMUMPS_INTR_STRUC` type).

### UNK-2 — C bridge struct layout — N/A under header-override approach

Original concern was that a hand-rolled `DMUMPS_STRUC_C` would have to
match the migrated Fortran derived type byte-for-byte. The header-
override approach (see B2 below) sidesteps this entirely: the C struct
is declared by upstream `dmumps_c.h` itself, and is INDEPENDENT from the
Fortran `DMUMPS_STRUC` (upstream's `dmumps_c()` decouples them via
field-by-field extraction). No layout match required.

### UNK-3 — `NO_SAVE_RESTORE` compile flag — RESOLVED 2026-05-02

`recipes/mumps.yaml` does not define `NO_SAVE_RESTORE` and the build
does not pass `-DNO_SAVE_RESTORE`, so the macro is **undefined** and
the `#if ! defined(NO_SAVE_RESTORE)` branches are active throughout
the staged tree (i.e. SAVE/RESTORE code is compiled in). Verified by
grepping `/tmp/stg-q/mumps/src/`:

```
$ grep -nE '#if[[:space:]]*!?[[:space:]]*defined\(NO_SAVE_RESTORE\)' \
       qmumps_save_restore.F qmumps_driver.F qini_defaults.F qmumps_f77.F
qmumps_save_restore.F:14:#if ! defined(NO_SAVE_RESTORE)
qmumps_f77.F:32: / :36: / :47: / :80: / :323:#if ! defined(NO_SAVE_RESTORE)
qmumps_driver.F:23 / :176 / :441 / :537 / :563 / :662 / :841 / :853 /
  :921 / :934 / :2065 — all `! defined(...)` form, all active
qini_defaults.F:682 / :1052 — same
mumps_print_defined.F:107: lone `defined(NO_SAVE_RESTORE)` (printed-config only)
```

Implication for the C bridge: the F77 prototype that the bridge calls
through (`${LIB_PREFIX}mumps_f77_`) is the SAVE/RESTORE-extended one,
matching upstream's default. No additional macro plumbing required.

### UNK-4 — Link-time satisfaction — RESOLVED 2026-05-02

Verified against the existing kind16 build under `/tmp/stg-q/build/`:

```
$ nm /tmp/stg-q/build/libqmumps-gfortran-13.a | grep -i pord | grep ' U '
(0 matches)
```

Zero unresolved PORD symbols in the migrated archive. The PORD bits
that MUMPS uses are reached through the C-side `mumps_pord.c` (compiled
in `mumps_c_bridge` for tests), not through Fortran callers. Linking
test executables with `${LIB_PREFIX}{mumps,scalapack,lapack,blas}` +
`MPI::MPI_Fortran` is sufficient — confirmed end-to-end by the kind16
suite already passing 23 mumps tests.

## Test design

### Approach

Differential precision against another MUMPS isn't feasible (no quad
reference exists). Instead:

1. Generate a small dense problem `A * x_true = b` at REAL(KIND=16) using
   `gen_*_quad` helpers (subset copied from `tests/lapack/common/`).
2. Convert to MUMPS triplet format (matrices small enough that all entries
   are kept; sparsity not required for entrypoint coverage).
3. Solve at the target precision with migrated MUMPS.
4. Cast the solution back to REAL(KIND=16) and check both:
   - **Forward error:** `tol_fwd  = 16  * n^3 * target_eps`  (matches
     `tests/lapack/linear_solve/test_dgesv.f90` n³ scaling — repo convention).
   - **Residual:**     `tol_resid = 64  * n   * target_eps`  (Higham, ASNA Thm 9.5,
     componentwise backward-error bound without growth factor).
5. Report each case via `prec_report` JSON; output lands in
   `<build>/precision_reports/<test_name>.<target>.json` (per-test-name to
   avoid clobbering under parallel CTest).

This satisfies the "always quad" testing convention — reference and all
comparisons happen at REAL(KIND=16), the migrated code runs at the target
precision.

### Coverage matrix — Fortran side

| File                       | Vary                                                                                              |
| -------------------------- | ------------------------------------------------------------------------------------------------- |
| `test_dmumps_basic`        | `n ∈ {8, 32}`, `SYM=0`, `JOB=6`, default ICNTL — smoke + tolerance baseline.                       |
| `test_dmumps_sym`          | `SYM ∈ {0, 1, 2}` with general / SPD / general-symmetric matrix construction.                      |
| `test_dmumps_jobs`         | Combined `JOB=6` vs phased `JOB=-1,1,2,3,-2`; compare **only** `INFOG(1)` and the solution `RHS` (other INFOG entries may legitimately differ). |
| `test_dmumps_nrhs`         | `NRHS ∈ {1, 3}`; per-rhs residual check; multi-rhs columns match single-rhs runs.                  |
| `test_dmumps_orderings`    | `ICNTL(7) ∈ {0(AMD), 2(AMF), 6(QAMD), 7(auto)}` — results equivalent within tolerance.             |
| `test_dmumps_errors`       | Invalid `N`, invalid `SYM`, mismatched `NNZ`, out-of-range entry → `INFOG(1) < 0`.                  |
| `test_zmumps_basic`        | `SYM=0`, `JOB=6` — smoke for complex.                                                              |
| `test_zmumps_sym`          | `SYM ∈ {0, 2}` (defer Hermitian PD `SYM=1` — needs special construction).                          |
| `test_zmumps_errors`       | Mirror of `dmumps_errors` for complex.                                                             |

ICNTL(1..4) set to `-1`/`0` to silence MUMPS output. Each test program calls
`MPI_INIT` / `MPI_FINALIZE` explicitly.

### C side — header-override bridge (replaces earlier hand-port plan)

Approach: ship a quad-precision shadow of upstream `mumps_c_types.h`
in `tests/mumps/c/include/`, put it FIRST on the include path, then
compile upstream `mumps_c.c`, `mumps_common.c`, `mumps_addr.c`, and
the upstream IO/save/thread/metis/scotch/pord helpers verbatim. See
B2 below for the resolved details and verified build/link recipe.

Why not the hand-port: with overrides, we ship ~80 lines of header
shadow versus ~800 lines of bridge code, and we automatically inherit
upstream's full C-side surface (assign/nullify/get pairs, OOC
helpers, save/restore, type-size utilities). When MUMPS upstream
upgrades, the override stays valid as long as the type-macro names
(`DMUMPS_REAL`, `ZMUMPS_COMPLEX`, ...) and the F77 symbol name
(`dmumps_f77_`) remain stable, which they have across MUMPS 5.x.

C drivers must call `MPI_Comm_c2f(MPI_COMM_WORLD)` (or use the
upstream `USE_COMM_WORLD` magic value `-987654`) for `comm_fortran`.
Upstream's `KEEP(500)=1` post-init logic and `NO_SAVE_RESTORE` flag
behavior are inherited unchanged.

### C-side coverage matrix (scoped, ~5 programs)

| File                       | Purpose                                                                       |
| -------------------------- | ----------------------------------------------------------------------------- |
| `test_dmumps_c_basic`      | `SYM=0`, `JOB=6` — smoke + residual via the bridge.                            |
| `test_dmumps_c_sym`        | `SYM ∈ {0,1,2}` — bridge handles all three.                                    |
| `test_dmumps_c_errors`     | Invalid input → `infog[0] < 0`.                                                |
| `test_dmumps_c_parity`     | Same problem solved through Fortran and C; outputs agree to `O(target_eps)`.   |
| `test_zmumps_c_basic`      | Complex smoke.                                                                |

## File layout (when implementation begins)

```
tests/mumps/
├── CMakeLists.txt               — gates: TARGET ${LIB_PREFIX}mumps,
│                                  HAVE_REAL16, reflapack_quad,
│                                  MPI_Fortran_FOUND + MPI_C_FOUND
├── README.md                    — scope + how-to-run (ctest via mpiexec)
├── TODO.md                      — this file
├── common/
│   ├── prec_kinds.f90           — copy from tests/lapack/common/
│   ├── compare.f90              — copy from tests/lapack/common/
│   ├── prec_report.f90          — adapt: filename = <test_name>.<target>.json
│   ├── test_data_mumps.f90      — sparse-from-dense quad helpers
│   └── target_mumps_body.fypp   — per-target wrapper template
├── target_kind16/
│   └── target_mumps.fypp        — kind16 wrapper module (only target wired today)
├── fortran/
│   ├── test_dmumps_basic.f90
│   ├── test_dmumps_sym.f90
│   ├── test_dmumps_jobs.f90
│   ├── test_dmumps_nrhs.f90
│   ├── test_dmumps_orderings.f90
│   ├── test_dmumps_errors.f90
│   ├── test_zmumps_basic.f90
│   ├── test_zmumps_sym.f90
│   └── test_zmumps_errors.f90
└── c/
    ├── include/
    │   ├── mumps_c_types.h      — quad-precision shadow (FIRST in -I)
    │   ├── mumps_int_def.h      — pin MUMPS_INTSIZE32
    │   ├── qmumps_c.h           — wraps dmumps_c.h with Q renames
    │   └── xmumps_c.h           — wraps zmumps_c.h with X renames
    ├── test_dmumps_c_basic.c
    ├── test_dmumps_c_sym.c
    ├── test_dmumps_c_errors.c
    ├── test_dmumps_c_parity.c
    └── test_zmumps_c_basic.c
```

## CMake gates (when implementation begins)

```cmake
if(NOT TARGET ${LIB_PREFIX}mumps)
    message(STATUS "tests/mumps: skipping — ${LIB_PREFIX}mumps not built")
    return()
endif()
if(NOT HAVE_REAL16)
    message(STATUS "tests/mumps: skipping — HAVE_REAL16 required")
    return()
endif()
if(NOT TARGET reflapack_quad)
    message(STATUS "tests/mumps: skipping — reflapack_quad not available")
    return()
endif()
if(NOT MPI_Fortran_FOUND OR NOT MPI_C_FOUND)
    message(STATUS "tests/mumps: skipping — MPI Fortran+C not found")
    return()
endif()
```

Test registration runs every executable via `mpiexec -n 1`:

```cmake
add_test(NAME mumps_${name}
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}>
                 ${MPIEXEC_POSTFLAGS}
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
set_tests_properties(mumps_${name} PROPERTIES RUN_SERIAL TRUE)
```

`RUN_SERIAL TRUE` is belt-and-braces for shared `precision_reports/` output;
drop after first runs confirm per-test-name JSON filenames don't collide.

## Carry-forward infrastructure blockers

### B5 — Migrator collapsed continuation-line inline `!` comments — RESOLVED 2026-04-29

When the migrator joined multi-line Fortran statements into a single
logical line for transform passes, inline `!` comments survived into the
joined string. After `reformat_fixed_line` re-split the joined string at
column 66, the `!` ended up mid-statement and swallowed every following
chunk. Original failure (qlr_core.F:1379):

```
+    0, ! L Panel          K, BLR_L)
```
→ `! L Panel` ate `K, BLR_L)`, breaking the call.

Fix in `pyengine/fortran_migrator.py`: `_strip_inline_comment()` helper
applied inside `_segment_fixed_form_statements` when joining
continuation lines. Inline comments are intentionally lost from the
joined / reformatted output (the price of correctness); single-line
statements with no transform keep their comments verbatim via the
per-line emit path.

Verified end-to-end:

```
$ cmake --build /tmp/stg-q/build --target qmumps -j8
…
[100%] Built target qmumps
```

`libmumps_common-gfortran-13.a` (~2 MB) and `libqmumps-gfortran-13.a`
(~22 MB) build cleanly. The same fix unblocks any other recipe whose
upstream Fortran has continuation-line inline comments.

### B6 — `FAC_FUTURE_NIV2_MOD` recipe gap — RESOLVED 2026-04-29

The integer-only module `MUMPS_FUTURE_NIV2` (defined in
`external/MUMPS_5.8.2/src/fac_future_niv2_mod.F`) was being classified
as PRECISION content and built into `${LIB_PREFIX}mumps`, but the
copy_files-resident `MUMPS_LOAD` (in `mumps_common`) `USE`s it.
mumps_common is built first → `mumps_load.F` failed with
`Cannot open module file 'mumps_future_niv2.mod'`.

Fixed by adding `FAC_FUTURE_NIV2_MOD` to `copy_files` in
`recipes/mumps.yaml`. (Note: `copy_files` matches uppercased filename
basenames, not module names — the entry is the file basename
`FAC_FUTURE_NIV2_MOD`, not the module name `MUMPS_FUTURE_NIV2`.)

### B7 — Two cmake/recipe trees diverge (fm-mumps vs fortran-migrator) — RESOLVED 2026-05-01

The split is gone. The mumps work that briefly lived in a sibling
`fm-mumps` tree was merged into fortran-migrator on the `tests`
branch (commits `9ec8357` "Merge branch 'mumps' into tests" and
`0b92df6` "Merge remote-tracking branch 'origin/mumps' into tests"),
and `fm-mumps` no longer exists as a separate working tree on this
machine. All three items the original entry called out as diverged
now sit in one place:

- `cmake/CMakeLists.txt` — section 9 (MUMPS) and section 10
  (libmpiseq) live in `cmake/CMakeLists.txt:747` and `:779`.
- `recipes/mumps.yaml` — `FAC_FUTURE_NIV2_MOD` is wired in
  `recipes/mumps.yaml:89`.
- `src/pyengine/__main__.py` — `('mumps', 'mumps.yaml')` is in
  `LIBRARY_ORDER` at line 724.

`pyengine stage` reads from this repo with no `--project-root`
override needed. The fm-mumps wrapper-script / recipe-drift CI
options listed in earlier drafts of this entry no longer apply —
keeping everything in one tree is the simpler outcome.

If a sibling MUMPS-only repo is ever introduced again, the right
move would still be option 2 (a wrapper script that pins
`--project-root` explicitly) over symlinks or a drift-check CI step.
But until then, this entry is closed.

### B2 — C interface (mumps_c.c, *mumps_c.h) not migrated — RESOLVED 2026-04-29 via header-override bridge

Recipe `skip_files` excludes every `*MUMPS_C` and `MUMPS_C_TYPES` header,
and `recipes/mumps.yaml` declares `language: fortran`. The C side is
provided in-tree by `tests/mumps/c/`:

- `tests/mumps/c/include/mumps_c_types.h` — quad-precision shadow of
  upstream's mumps_c_types.h. Sets `DMUMPS_REAL`/`ZMUMPS_REAL` to
  `__float128` and `mumps_double_complex` to a struct of two
  `__float128`. When the test build's include path lists
  `tests/mumps/c/include/` BEFORE `external/MUMPS_5.8.2/include/`,
  every consumer of `<mumps_c_types.h>` (upstream `mumps_c.c`,
  `mumps_common.c`, `dmumps_c.h`, `zmumps_c.h`, ...) sees the quad
  typedefs.
- `tests/mumps/c/include/mumps_int_def.h` — pin to `MUMPS_INTSIZE32`
  (upstream's build script generates this from a template; we don't
  run that script here).
- `tests/mumps/c/include/qmumps_c.h` and `xmumps_c.h` — small wrappers
  around upstream `dmumps_c.h` / `zmumps_c.h` that macro-rename the
  user-visible entry points (`dmumps_c → qmumps_c`,
  `DMUMPS_STRUC_C → QMUMPS_STRUC_C`, `dmumps_f77_ → qmumps_f77_` and
  the Z analogues).

Build glue (to add when `tests/mumps/CMakeLists.txt` lands):

```cmake
set(_MUMPS_UPSTREAM ${PROJECT_SOURCE_DIR}/external/MUMPS_5.8.2)

# 1. Type-agnostic C runtime (upstream sources, compiled verbatim).
add_library(mumps_c_runtime STATIC
    ${_MUMPS_UPSTREAM}/src/mumps_common.c
    ${_MUMPS_UPSTREAM}/src/mumps_addr.c
    ${_MUMPS_UPSTREAM}/src/mumps_io.c
    ${_MUMPS_UPSTREAM}/src/mumps_io_basic.c
    ${_MUMPS_UPSTREAM}/src/mumps_io_err.c
    ${_MUMPS_UPSTREAM}/src/mumps_io_thread.c
    ${_MUMPS_UPSTREAM}/src/mumps_save_restore_C.c
    ${_MUMPS_UPSTREAM}/src/mumps_register_thread.c
    ${_MUMPS_UPSTREAM}/src/mumps_thread.c
    ${_MUMPS_UPSTREAM}/src/mumps_thread_affinity.c
    ${_MUMPS_UPSTREAM}/src/mumps_numa.c
    ${_MUMPS_UPSTREAM}/src/mumps_flytes.c
    ${_MUMPS_UPSTREAM}/src/mumps_config_file_C.c
    ${_MUMPS_UPSTREAM}/src/mumps_metis.c
    ${_MUMPS_UPSTREAM}/src/mumps_metis64.c
    ${_MUMPS_UPSTREAM}/src/mumps_metis_int.c
    ${_MUMPS_UPSTREAM}/src/mumps_pord.c
    ${_MUMPS_UPSTREAM}/src/mumps_scotch.c
    ${_MUMPS_UPSTREAM}/src/mumps_scotch64.c
    ${_MUMPS_UPSTREAM}/src/mumps_scotch_int.c
)
target_include_directories(mumps_c_runtime PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/c/include
    ${_MUMPS_UPSTREAM}/include
    ${_MUMPS_UPSTREAM}/src)
target_compile_definitions(mumps_c_runtime PRIVATE Add_)
target_link_libraries(mumps_c_runtime PUBLIC MPI::MPI_C)

# 2. Per-arithmetic dispatch — compile mumps_c.c twice with overrides.
add_library(qmumps_c_bridge OBJECT ${_MUMPS_UPSTREAM}/src/mumps_c.c)
target_include_directories(qmumps_c_bridge PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/c/include
    ${_MUMPS_UPSTREAM}/include ${_MUMPS_UPSTREAM}/src)
target_compile_definitions(qmumps_c_bridge PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_d Add_
    dmumps_f77_=qmumps_f77_
    dmumps_set_tmp_ptr_=qmumps_set_tmp_ptr_
    dmumps_c=qmumps_c
    DMUMPS_STRUC_C=QMUMPS_STRUC_C)

add_library(xmumps_c_bridge OBJECT ${_MUMPS_UPSTREAM}/src/mumps_c.c)
target_include_directories(xmumps_c_bridge PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/c/include
    ${_MUMPS_UPSTREAM}/include ${_MUMPS_UPSTREAM}/src)
target_compile_definitions(xmumps_c_bridge PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_z Add_
    zmumps_f77_=xmumps_f77_
    zmumps_set_tmp_ptr_=xmumps_set_tmp_ptr_
    zmumps_c=xmumps_c
    ZMUMPS_STRUC_C=XMUMPS_STRUC_C)

add_library(mumps_c_bridge STATIC
    $<TARGET_OBJECTS:qmumps_c_bridge>
    $<TARGET_OBJECTS:xmumps_c_bridge>)
target_link_libraries(mumps_c_bridge PUBLIC mumps_c_runtime)
```

Note on the renames:

- `dmumps_f77_ → qmumps_f77_` — Fortran symbol called *from* C; rename
  so the bridge resolves to the migrated archive's entry.
- `dmumps_set_tmp_ptr_ → qmumps_set_tmp_ptr_` — same: another Fortran
  symbol called from C.
- `dmumps_set_tmp_ptr_c_` is **NOT** renamed — it's a C symbol called
  *from* the migrated Fortran, which still calls it under the original
  D name (per the migrator's keep-kind manifest convention).
- Same pattern: every `dmumps_assign_*` / `dmumps_nullify_c_*` defined
  by `mumps_c.c` keeps its D prefix because the migrated qmumps_f77_
  Fortran calls it that way.

Verified end-to-end on 2026-04-29 by hand-linking a JOB=-1/-2
roundtrip:

```
$ mpiexec -n 1 ./test_qmumps_init
Entering QMUMPS 5.8.2 from C interface with JOB =  -2
      executing #MPI =      1, without OMP
```

Distributed mapping / permutation retrieval (sym_perm, uns_perm, etc.)
work the same way they do upstream — our bridge just re-uses upstream
`mumps_common.c` which provides the assign/nullify/get pairs.

### B3 — Sequential MPI stub (libmpiseq) — PARTIALLY RESOLVED 2026-04-29

`libmpiseq.a` is now buildable. `external/MUMPS_5.8.2/libseq/` is
copied to `_mpiseq_src/` by `pyengine stage` (added to `_std_dirs` in
`src/pyengine/__main__.py`), and `cmake/CMakeLists.txt` (section 10)
declares `add_library(mpiseq STATIC mpi.f mpic.c elapse.c)` with the
required `-DAdd_` Fortran-symbol-mangling flag.

**Build verified:**

```
$ cmake --build /tmp/stg-q/build --target mpiseq
[..]
[100%] Linking Fortran static library libmpiseq.a
[100%] Built target mpiseq
```

`libmpiseq.a` exports 157 symbols — the full MPI / BLACS / ScaLAPACK
stub surface upstream MUMPS uses for serial operation.

**Supplementary C-MPI stubs:** `tests/mumps/c/mpiseq_c_stubs.c`
provides ~30 additional C-side `MPI_*` functions (Send, Recv, Isend,
Bcast, Type_vector, ...) referenced by the standard BLACS / PBLAS C
archives. Compiled separately and linked alongside libmpiseq.

**Open: full sequential link is non-trivial.** The migrated qmumps
archive depends on Q-prefixed Fortran routines (`pqgetrs_`,
`pqgetrf_`, `pqpotrf_`, `pqpotrs_`) that libmpiseq provides only as
the upstream D-prefixed equivalents (`pdgetrs_`, ...). And keeping
`libqscalapack`/`libqblacs`/`libptzblas` in the link line for those
Q-symbols pulls in duplicate definitions of the precision-agnostic
descriptor / utility routines (`descinit_`, `numroc_`, `infog2l_`,
`indxg2p_`, `chk1mat_`, `pchk2mat_`, `pxerbla_`, `descset_`,
`blacs_gridinit_`, `blacs_gridinfo_`, `blacs_gridexit_`) — those are
defined by BOTH libmpiseq's `mpi.f` and the standard scalapack/blacs
archives.

Fix would require either (a) building a stripped variant of
`libqscalapack` / `libqblacs` that omits the routines libmpiseq
already covers, or (b) extending libmpiseq's `mpi.f` with Q/X
analogues of `pdgetrs_` etc. Both are sizable work. Tests today
keep using `mpiexec -n 1` against real MPI (verified end-to-end).
libmpiseq is left available as a future option for environments
without MPI installed.

**Update 2026-05-02 (path b, partial):** Q/X/E/Y per-precision stubs
landed in `cmake/mpiseq_qx_stubs.f` and are appended to the `mpiseq`
target via `cmake/CMakeLists.txt:835`. Twenty new stubs exported
(`p[qxey]{getrf,getrs,potrf,potrs,trtrs}_`), each cloning upstream
`mpi.f`'s "should not be called / STOP" body. Authored in-tree
because `external/` is read-only.

Multifloats (M/W) prefixes intentionally NOT covered yet — adding
sequential stubs without the matching test-harness work (Group A /
B9b) would be premature. Adding `pmgetrf_` etc. is a follow-up to
B9b whenever multifloats Fortran tests come online.

The duplicate-symbol concern from path (b) does not bite here:
upstream `mpi.f` defines only PD* — the Q/X/E/Y names are fresh,
so adding them does not collide with `libqscalapack` /
`libqblacs` either, and a stripped-variant recipe is no longer the
prerequisite that path (a) implied.

Linking a test executable against `mpiseq` (instead of
`MPI::MPI_Fortran`) plus `libqscalapack`/`libqblacs` should now
resolve all symbols for kind10 + kind16, but no test in this repo
exercises that path today (mpiexec -n 1 stays the default). Spot-
check via `nm libmpiseq.a | grep ' T pqgetrf'`.

### B4 — Only `kind16` target supported  *(PARTIAL — kind10 wired 2026-04-30)*

**Progress 2026-04-30**: kind10 EP overrides authored
(`recipes/mumps/kind10/{mumps_memory_mod_ep.F,mumps_lr_stats_ep.F}`),
recipe wired, `emumps-gfortran-13.a` builds clean against staged
mumps tree. `tests/mumps/target_kind10/target_mumps.fypp` added.
A latent line-length issue surfaced (kind10's
`MPI_C_LONG_DOUBLE_COMPLEX` rewrite pushes fixed-form `.F` past
column 72) and was fixed by adding `-ffixed-line-length-none`
(gfortran/flang) / `-extend-source` (Intel) to
`add_migrated_fortran_library` in `cmake/CMakeLists.txt`.

Multifloats overrides also authored
(`recipes/mumps/multifloats/{...}.F`) and recipe-wired but the
build hits a separate, pre-existing migrator bug (see B8 below).

Remaining for B4: multifloats build (blocked on B8); `INFOG(20)`
sanity-check across targets — multifloats is now testable as of
2026-05-02 (B9b RESOLVED), so the cross-target sanity check can run
on all three.

**Update 2026-05-02:** `tests/mumps/fortran/test_dmumps_infog20.f90`
landed for kind10 + kind16. The test factors a fixed n=32 SYM=0
JOB=6 problem, reads `id%INFOG(20)`, and asserts the result lies in
the structural bound `[n², 50·n²]`. The structural bound is wider
than the ±5% kind16-baseline check originally specified — it ships
before any per-target baseline has been captured, and catches the
gross-mis-sizing failure mode (off-by-2× / off-by-4× from a wrong
override constant) without needing a precaptured number. Tighten
to ±5% once two builds have logged the actual `INFOG(20)` values
into `precision_reports/test_dmumps_infog20.<target>.json` and the
spread is known.

Multifloats coverage on this test is deferred: ctest skips test
executables for `multifloats`
(`tests/mumps/CMakeLists.txt:52-59`'s `_MUMPS_TESTS_SKIP_EXECUTABLES`),
and Group A (B9b) is the path that re-enables them. Re-evaluate
once Group A lands and pull multifloats into the cross-target ±5%
check at that point.

Recipe `overrides:` in `recipes/mumps.yaml` declares hand-written
extended-precision substitutes only for the kind16 target:

```yaml
overrides:
  kind16:
    src_dir: mumps
    files:
      - mumps_memory_mod_ep.F
      - mumps_lr_stats_ep.F
```

The other targets in this repo (`kind10`, `multifloats`) have no
analogous files, so the recipe falls back to the upstream DP source
for those modules — which encodes byte-size and stride constants
based on `sizeof(REAL(KIND=8))` and is therefore wrong at any other
precision. End result: `pyengine stage --target kind10 --libraries mumps`
might link but produces a binary with mis-sized internal records;
`multifloats` likely fails outright.

#### Why these two specific files need overrides

Upstream MUMPS uses `MUMPS_MEMORY_MOD` and `MUMPS_LR_STATS` to
record peak memory and low-rank statistics in fixed-width records.
Several constants compute byte sizes from the working precision via
hardcoded literals like `8` (DP real) or `16` (DP complex). The
migrator promotes `REAL(KIND=8)` declarations to the target kind
correctly, but it can't translate "8 bytes per real" into "16 bytes
per quad" — that's a numeric literal in arithmetic context, not a
kind parameter. Hence the hand-written EP file with pre-computed
literals: upstream's `8` → override's `16`, etc.

Other modules in MUMPS (`MUMPS_LOAD`, `MUMPS_OOC_*`, etc.) get
through the migrator unscathed because they either avoid
`sizeof()`-style arithmetic or use `STORAGE_SIZE` / kind queries
that auto-adjust.

#### What "fix B4" looks like

1. **Author per-target override files** in `recipes/mumps/`. Easiest
   path: copy the existing `_ep.F` files as a skeleton, then walk
   through every numeric literal and confirm/adjust against the
   actual byte-width of the target real type:

   - **kind10** (x87 `long double`, 16 bytes on x86-64 with 16-byte
     alignment per Intel ABI; complex = 32 bytes):
     `mumps_memory_mod_kind10.F`, `mumps_lr_stats_kind10.F`
   - **multifloats kind16** (double-double, 16 bytes per real, 32
     per complex — same byte sizes as kind16, but the SIMD path may
     pack differently; compare carefully against the existing
     `_ep.F` to spot divergence):
     `mumps_memory_mod_multifloats.F`, `mumps_lr_stats_multifloats.F`
   - **multifloats kind32** (quad-double, 32 bytes per real, 64
     per complex) — only relevant if multifloats kind32 becomes a
     supported target.

2. **Register them** under per-target sections in
   `recipes/mumps.yaml`. Sibling recipes (`recipes/blacs.yaml`,
   `recipes/pblas.yaml`, `recipes/scalapack_c.yaml`) already use
   per-target overrides for `multifloats:` with `src_dir:
   <lib>/mfc_overrides`; follow that pattern. Result:

   ```yaml
   overrides:
     kind16:
       src_dir: mumps
       files: [mumps_memory_mod_ep.F, mumps_lr_stats_ep.F]
     kind10:
       src_dir: mumps
       files: [mumps_memory_mod_kind10.F, mumps_lr_stats_kind10.F]
     multifloats:
       src_dir: mumps/mfc_overrides
       files: [mumps_memory_mod_multifloats.F, mumps_lr_stats_multifloats.F]
   ```

3. **Build smoke per target.** After staging, rebuild `qmumps` →
   `emumps` (kind10) → `wmumps` (multifloats kind16) — the migrator
   prefix-replaces D with the target letter (E/Y for kind10, M/W for
   multifloats per `project_migrator_prefix_convention`).

4. **Sanity-check `INFOG(20)`** ("number of real-precision words
   used") and the peak-memory report after factorization. Those are
   the most direct user-visible signals that the byte-width math is
   correct. If `INFOG(20)` is wildly off versus an equivalent kind16
   run scaled by precision ratio, the override constants are wrong.

5. **Drop in test wrappers** (G2 below) — once the build works, the
   tests/mumps/ directory needs `target_kind10/target_mumps.fypp`
   and `target_multifloats/target_mumps.fypp` mirroring the existing
   kind16 wrapper.

#### Existing kind16 override is also incomplete

The current `recipes/mumps/mumps_memory_mod_ep.F` does NOT mirror
the upstream module's structure for size constants. Upstream
`mumps_memory_mod.F:23` declares module-level

```fortran
INTEGER(8), PRIVATE :: ISIZE, I8SIZE, SSIZE, DSIZE, CSIZE, ZSIZE
```

populated at startup by `MUMPS_MEMORY_SET_DATA_SIZES` via
`MUMPS_SIZE_C(D(1), D(2), DSIZE)` etc., and then used everywhere
throughout the module as `int(MINSIZE,8)*DSIZE` for byte-accounting.

The override sidesteps this by using a per-subroutine local
`INTEGER(8) :: ELSIZE` populated as
`ELSIZE = int(storage_size(1.0_16),8)/8_8` inside each of
`MUMPS_QREALLOC` / `MUMPS_XREALLOC`. Functional for the override's
own bookkeeping, but it means there are **no `QSIZE` / `XSIZE`
module-level constants** the way upstream has `DSIZE` / `ZSIZE`,
and the override is therefore not a drop-in structural mirror.

The future overrides need to mirror the standard-precision module
faithfully:

- kind10 override needs **`ESIZE`** (kind10 real) and **`YSIZE`**
  (kind10 complex) — module-level `INTEGER(8), PRIVATE`,
  populated by an extended `MUMPS_MEMORY_SET_DATA_SIZES`.
- kind16 override needs **`QSIZE`** and **`XSIZE`** added —
  retroactive fix to bring the existing `_ep.F` into structural
  parity. (Today it works because the local `ELSIZE` covers
  internal byte accounting, but a module-level `QSIZE` is the
  pattern.)
- multifloats override needs **`MSIZE`** (multifloats real) and
  **`WSIZE`** (multifloats complex).

See the `project_migrator_prefix_convention` memory for the
precision-letter mapping (E/Y for kind10, Q/X for kind16, M/W for
multifloats).

#### Why deferred

No active consumer is waiting for kind10 / multifloats MUMPS, and
sizing the records correctly requires careful inspection of the
upstream `mumps_memory_mod.F` to enumerate every precision-dependent
literal — the override is a one-shot hand-port and getting it wrong
is silent (mis-sized buffers in statistics records, not factorization
errors). The other libraries in this repo (BLAS / LAPACK /
ScaLAPACK / xBLAS) work across all three target families today;
MUMPS is the lone holdout until someone needs a non-quad sparse
solver.

#### Files to read when picking this up

- `recipes/mumps/mumps_memory_mod_ep.F` — kind16 reference; diff
  vs upstream `external/MUMPS_5.8.2/src/mumps_memory_mod.F` to see
  exactly which literals were changed.
- `recipes/mumps/mumps_lr_stats_ep.F` — kind16 reference for LR
  statistics module.
- `external/MUMPS_5.8.2/src/mumps_memory_mod.F` — upstream DP source.
- `external/MUMPS_5.8.2/src/mumps_lr_stats.F` — upstream LR source.
- `recipes/blacs.yaml`, `recipes/pblas.yaml` — examples of recipes
  with per-target overrides for multifloats.

## Open coverage gaps

### G1 — Fortran/C cross-language parity test — RESOLVED 2026-04-30

Landed as `fortran/test_dmumps_c_parity.f90` and `test_zmumps_c_parity.f90`
plus `c/c_parity_helpers.c` (the Fortran-callable bridge into
qmumps_c / xmumps_c). Each test drives the same generated problem
through both paths and verifies bit-identical solutions. Catches
silent struct-extraction corruption that the per-side basic tests
can't (since they only check `result ≈ x_true` to ~33 digits, not
that the two paths agree exactly).

Avoided the `__float128` ISO_C_BINDING gymnastics by using
`real(16)` / `complex(16)` in the Fortran interface block — gfortran
and gcc agree on the by-reference layout for those kinds, so the
bridge function declared as `__float128 *` / `mumps_double_complex *`
on the C side links cleanly.

Verified `fortran-vs-c = 0.0e+00` on both D and Z paths.

### G2 — `target_kind10/` / `target_multifloats/` placeholder dirs — RESOLVED 2026-05-02

Both wrapper directories now exist on disk:

- `tests/mumps/target_kind10/target_mumps.fypp`
- `tests/mumps/target_multifloats/target_mumps.fypp`

Each ships a precision-aware `c/include/mumps_c_types.h` shadow plus
the per-target `*mumps_c.h` wrappers (B9 cleanup). Closed.

## Defects discovered during runs

### D1 — MUMPS 5.8.2 doesn't validate out-of-range matrix entries / negative N

Discovered 2026-04-30 while writing `test_dmumps_errors.f90` /
`test_zmumps_errors.f90`. Inputs that the test plan expected MUMPS to
reject with `INFOG(1) < 0` instead crash inside the analysis or
factorization phase:

- `id%N = -1` → SIGSEGV before MUMPS prints any diagnostic.
- `id%IRN(1) = N + 5` or `id%JCN(1) = N + 3` → SIGSEGV inside
  `qmumps_validate_input_` or its callees.
- Mismatched `NNZ` (e.g. `NNZ = 0` with non-empty IRN/JCN/A) → SIGSEGV
  during matrix reformatting.

The MUMPS 5.8.2 user manual documents `INFOG(1) = -16` (N out of range)
and `INFOG(1) = -6` (structurally singular) as the intended behavior.
The migrated qmumps inherits whatever validation upstream actually has,
which doesn't cover these cases. The test files for error coverage
were therefore dropped — keeping them would test MUMPS-side robustness
that doesn't exist, not migrator correctness. Re-introduce when MUMPS
upstream tightens validation, or write a lightweight wrapper that
sanitizes inputs before calling `qmumps_c`.

**Wrapper landed 2026-05-02.** `target_mumps` now exports
`check_dmumps_input` / `check_zmumps_input` plus `MIC_*` codes
(`tests/mumps/common/target_mumps_body.fypp`). Restored tests live at
`tests/mumps/fortran/test_dmumps_errors.f90` and
`tests/mumps/fortran/test_zmumps_errors.f90`; they exercise BAD_N
(neg + zero), BAD_NNZ, BAD_IRN (high + zero), BAD_JCN (high + zero),
SIZE_MISMATCH, plus a final valid-input pass that reaches MIC_OK and
factors via JOB=6. The wrapper layers above MUMPS — does not touch
`id%INFOG(1)`. Tests target-portable via `target_qmumps`/`target_xmumps`.

## New follow-ups (2026-04-30, surfaced during B4 partial completion)

### B8 — Migrator's `USE multifloats` injector mishandles CPP-interspersed argument continuation — RESOLVED 2026-05-01

Resolved by tracking parenthesis depth across the SUBROUTINE/FUNCTION
header walker in `insert_use_multifloats` (commit `bf33cde`). The walker
also now treats CPP `#if/#endif/#else/#elif/#define` directives as
transparent. Helper `_count_open_parens` returns the net paren delta
for a line, ignoring quoted strings and inline comments.

Two follow-up multifloats mmumps build issues surfaced after B8 was
fixed; these are tracked as B8b and B8c below.

#### Original B8 description (for reference)

Surfaced 2026-04-30 attempting to build `mmumps`. Migrated files like
`/tmp/stg-mf-mumps/mumps/src/mana_aux.F:26-38` show the multifloats
USE statement injected INSIDE the formal argument list, between the
last regular continuation line and the post-CPP continuation:

```fortran
      SUBROUTINE MMUMPS_ANA_F(N, NZ8, IRN, ICN, LIWALLOC, ...
     +CNTL4, COLSCA, ROWSCA
      USE multifloats, only: ...        <-- INJECTED HERE (wrong)
#if defined(metis) ...
     &          , METIS_OPTIONS
#endif
     &          , NORIG_ARG, ..., GCOMP
     & )
      USE MUMPS_ANA_ORD_WRAPPERS         <-- legitimate USE site
      ...
```

The injector treats the first non-continued line ending in something
other than `)` as the SUBROUTINE statement terminator and inserts USE
right after, missing that the formal arg list resumes after a CPP
`#if/#endif` block via `&` continuation and only closes at `& )`
several lines down. The injection ends up inside the formal arg list,
breaking compilation.

Affects: every multifloats `.F` file that has CPP-mediated continuation
inside a SUBROUTINE/FUNCTION header. Many MUMPS files trip this because
they conditionalize METIS/parmetis arguments inline.

Workaround: none today — fully blocks `mmumps` build. Engine fix
should treat CPP `#if/#endif` blocks as transparent for the
"end-of-statement" search inside SUBROUTINE/FUNCTION headers, or use
a Fortran-aware parser instead of textual scanning at the injection
site (`fortran_migrator.py` — multifloats post-pass).

### B8b — `int(real64x2_expr, KIND_INT)` has no multifloats overload — RESOLVED 2026-05-01

Surfaced 2026-05-01 once B8 unblocked the SUBROUTINE-header parse for
`mana_aux.F`. The migrator's `replace_intrinsic_calls` wraps the inner
of `INT(double_expr, 8)` calls with `real64x2(...)` (so DOUBLE PRECISION
arguments lift to multifloats), producing `int(real64x2(...), 8)`. But
multifloats's `int` interface is only `dd_int(real64x2) -> integer`
(default kind) — there is no `int(real64x2, kind=8)` overload, and
gfortran refuses with "Generic function 'int' is not consistent with a
specific intrinsic interface" (mana_reordertree.F:623,
mumps_ooc.F:215/219/224/228 etc.).

Resolved by `_rewrite_int_kind_on_real64x2` in `fortran_migrator.py`:
when the arg of `INT(...)` contains a `real64x2(` token AND the call
has a kind specifier (numeric or `KIND=N`), wrap the arg with
`dd_to_double(...)` so the standard Fortran INT intrinsic handles the
kind selector. `dd_to_double` is a public multifloats helper that
returns `real(dp)`. Values that flow through this pattern (cost
counters, matrix size estimates) fit comfortably in double precision.

### B8c — `MODULE PROCEDURE` line falsely matched as procedure header — RESOLVED 2026-05-01

Surfaced 2026-05-01 building `mmumps_ana_aux_par_m` after B8b. In
`mana_aux_par.F:25-31`:

```fortran
INTERFACE MMUMPS_ANA_F_PAR
MODULE PROCEDURE MMUMPS_ANA_F_PAR
USE multifloats, only: ...        <-- INJECTED HERE (wrong)
END INTERFACE
```

The `_PROC_HEADER_RE` in `fortran_migrator.py` matched the bare keyword
``MODULE`` with `\b` boundary, so ``MODULE PROCEDURE`` lines inside
INTERFACE blocks were treated as procedure-header starts and had the
USE clause inserted right after — but USE inside an INTERFACE block
between MODULE PROCEDURE and END INTERFACE is illegal Fortran.

Resolved by adding a negative lookahead `(?!\s+PROCEDURE\b)` after
``MODULE`` in both `_PROC_HEADER_RE` (line 1832) and the local
`proc_header_re` inside `insert_use_multifloats` (line 2030).

### B8e — `MPI_FLOAT64X2` undeclared in migrated MUMPS Fortran — RESOLVED 2026-05-01

`external/multifloats-mpi/multifloats_mpi_f.f90` ships a Fortran module
that bind-to-C declares each handle as a default-INTEGER variable, and
`external/multifloats-mpi/multifloats_mpi.cpp` populates them via
`MPI_Type_c2f` / `MPI_Op_c2f` inside `multifloats_mpi_init`. The
migrator's new `insert_use_multifloats_mpi_f` post-pass walks each
procedure header (with proper depth tracking through INTERFACE-block
inner SUBROUTINEs) and emits a `USE multifloats_mpi_f` line whenever
the body references one of the handle names.

Scope: real datatype `MPI_FLOAT64X2` + complex datatype
`MPI_COMPLEX64X2` (renamed from the original `MPI_COMPLEX128X2` to
match the upstream `cmplx64x2` type) + all reduction ops
(`MPI_DD_*`, `MPI_ZZ_*`) are exposed. The complex datatype handle was
added in the B8g cleanup commit when wmumps multifloats build was
brought up alongside mmumps.

After this fix all multifloats mumps source files (m-prefix and
w-prefix) compile end-to-end through `_rewrite_mpi_datatypes` +
`USE multifloats_mpi_f` + the bind-to-C handles. A handful of unrelated
migrator bugs that the previously-blocking files masked then surfaced
— see B8f and B8h below.

### B8f — Pre-existing migrator pipeline bugs surfaced once B8e unblocked mmumps

Once B8e let the build progress through the multifloats mumps tree, a
cluster of latent migrator bugs became visible. Resolved in the same
commit as B8e:

- `insert_use_multifloats`'s decl-block walker did not list `INCLUDE`
  among declaration keywords, so PARAMETER-derived assignments landed
  between `IMPLICIT NONE` and `INCLUDE 'mpif.h'`. Same walker also
  stopped at keep-kind sentinels (`__KEEPKIND_DP__` etc.) since they
  no longer matched the bare `DOUBLE PRECISION` keyword. Both added
  to the whitelist.
- `convert_parameter_stmts` lost the `#if defined(...)` guard when
  converting a guarded PARAMETER to a runtime assignment, so files
  like `wfac_mem_stack_aux.F` ended up with `ZERO = ...` outside the
  `#if defined(ZERO_TRIANGLE)` block that declared `ZERO`. Now wraps
  the assignment in the same guard.
- `_rewrite_int_kind_on_real64x2` emitted `dd_to_double(...)` and
  `_rewrite_int_of_complex` emitted `DD_REAL(...)`. Neither symbol
  exists in upstream multifloats — the multifloats public generics
  are `dble` and `real`. Fixed to emit those; `_build_use_only_clause`
  predicts the rewrite (when `int(` and `real64x2(` co-occur, add
  `dble` to referenced) so the only-clause is in place.
- The depth-tracking `end_proc_re` used by
  `insert_use_multifloats_mpi_f` falsely matched `END IF` / `END DO`
  (the keyword group was optional and `\w*` swallowed anything).
  Now requires a procedure keyword whitelist after `END` if any
  word follows.

### B8g — Migrator wraps PARAMETER RHS literal even when LHS type is keep-kind DOUBLE PRECISION — RESOLVED 2026-05-01

`replace_literals` now early-returns the line unchanged in
``constructor`` mode whenever any keep-kind sentinel
(``__KEEPKIND_DP__`` / ``__KEEPKIND_DBLE__`` / ``__KEEPKIND_DCMPLX__``)
appears on the line. The keep-kind LHS / wrapper preserves the
original DP semantics, so the RHS DP literal must stay un-wrapped to
avoid the ``Cannot convert TYPE(real64x2) to REAL(8)`` mismatch. The
``kind_suffix`` path is unaffected (REAL(KIND=N) → REAL(8) implicit
narrowing is legal).

### B8h — Pre-existing migrator bugs surfaced once B8g unblocked mmumps multifloats — RESOLVED 2026-05-01

Resolved in the same commit as B8g. The mmumps Fortran build for the
multifloats target now compiles end-to-end (real m-prefix + complex
w-prefix). Items in this cluster:

- **Combined-form ``TYPE, PARAMETER ::`` declarations** were left as
  invalid PARAMETER initializers (``cmplx64x2(...)`` etc. is not a
  constant expression). ``convert_parameter_stmts`` now matches the
  combined form (``COMPLEX(kind=8), PARAMETER :: ZERO = (0.0,0.0)``,
  ``DOUBLE PRECISION, PARAMETER :: DZERO = 0.0D0``, ``TYPE(...)``,
  ``INTEGER(...)``) and splits it into a plain declaration + runtime
  assignment. The INTEGER + complex-literal degenerate form (upstream
  ``INTEGER, PARAMETER :: ZERO = (0.0D0,0.0D0)`` in zsol_fwd_aux.F:1095)
  is detected and the type is promoted to the multifloats complex type.
- **File-scope COMPLEX/REAL ambiguity for known-constant names**.
  ``strip_known_constants_from_decls`` and
  ``_filter_known_constants_from_decl`` now accept a per-file
  ``complex_names`` set and skip stripping a known-constant name that
  is also declared as COMPLEX in some other procedure of the same
  file. Without this, a global ``ZERO -> DD_ZERO`` rename mistyped
  the COMPLEX scope.
- **Renames overwriting LHS of converted-PARAMETER assignments**.
  ``replace_known_constants`` now detects ``^name = ...`` lines whose
  LHS would be renamed by ``removed_known``; the LHS rename is
  suppressed for those lines so the runtime assignment stays
  ``ZERO = ...`` rather than becoming ``DD_ZERO = ...`` (an
  assignment to a multifloats public PARAMETER, which gfortran rejects
  with ``Named constant in variable definition context``).
- **Keep-kind sentinel masking declaration detection.**
  ``_DECL_KEYWORD_RE`` now matches ``__KEEPKIND_DP__`` so
  ``replace_known_constants`` correctly skips declaration lines
  whose ``DOUBLE PRECISION`` token has been masked.
- **``int(var, KIND_INT)`` over real64x2 variables.**
  ``_rewrite_int_kind_on_real64x2`` previously triggered only on the
  literal-wrapper form ``int(real64x2(...), 8)``; it now also fires
  when the inner expression is a single identifier that the file scan
  recognises as ``TYPE(real64x2)``. Recognises ``int(...)`` with a
  preceding space (``int (`` in continuation-merged statements).
- **INTERFACE blocks treated as procedure scopes.**
  ``insert_use_multifloats``, ``_scope_index_at`` and
  ``specialize_use_module`` now track INTERFACE depth: prototype
  SUBROUTINE/FUNCTION declarations inside an INTERFACE block are NOT
  counted as new scopes for runtime-assignment placement, but the
  decl-block walker injects ``USE <module>`` after each inner header
  so the prototype body sees ``real64x2`` / ``cmplx64x2`` (interface
  scopes do not inherit the enclosing procedure's USE clauses).
  Runtime assignments are placed AFTER the entire decl section
  (including any INTERFACE blocks), not before it.

### B8d — `dble(...)` masked by keep-kind sentinel missing from USE only-list — RESOLVED 2026-05-01

Surfaced 2026-05-01 after B8b. `recipes/mumps/keep-kind.manifest`
includes `dana_aux.F:4291` to preserve a `dble(PEAK)` call to a
verbatim (copy_files) callee `MUMPS_DISTRIBUTE` that takes a real(8).
The keep-kind sentinel masks `dble(` as `__KEEPKIND_DBLE__(` before
migration, so `_scan_referenced_identifiers` never sees `dble` and
`_build_use_only_clause` omits it from the only-list. After the
sentinel restore swaps `dble` back, the call site dispatches to
gfortran's intrinsic `dble`, which doesn't accept the multifloats
`real64x2` type ("'a' argument of 'dble' intrinsic at (1) must have
a numeric type").

Resolved in `_build_use_only_clause`: detect `_KK_DBLE_SENTINEL` in the
proc body text and unconditionally add `dble` to the referenced set so
multifloats's `dd_dble(real64x2) -> real(dp)` overload gets imported.

### B9 — C-bridge / test infrastructure hardcoded to kind16 — RESOLVED 2026-05-01

The bridge cmake glue (tests/mumps/CMakeLists.txt) and the per-
arithmetic header shadow (tests/mumps/target_${TARGET_NAME}/c/include/)
now follow `${LIB_PREFIX}` / `${LIB_PREFIX_COMPLEX}`. Each target
ships its own:

- `mumps_c_types.h` — precision-aware DMUMPS_REAL / ZMUMPS_COMPLEX
  typedefs (kind16: `__float128` / `mumps_double_complex` of two
  `__float128`; kind10: `long double` / two `long double`;
  multifloats: `mumps_float64x2` POD / two such with `.r/.i`).
- `${LIB_PREFIX}mumps_c.h` / `${LIB_PREFIX_COMPLEX}mumps_c.h` —
  thin wrappers around upstream `dmumps_c.h` / `zmumps_c.h` with
  the per-target macro renames.

The cmake bridge OBJECT-library names are now
`${LIB_PREFIX}mumps_c_bridge` / `${LIB_PREFIX_COMPLEX}mumps_c_bridge`
(replacing the old fixed `qmumps_c_bridge`/`xmumps_c_bridge`), with
`-Ddmumps_c=${LIB_PREFIX}mumps_c` etc. driving the per-arithmetic
dispatch. `c_parity_helpers.c` is templated through the C macros
`TARGET_REAL_HEADER` / `TARGET_REAL_MUMPS_C` / `TARGET_REAL_STRUC_C`
(plus complex counterparts) so it builds unchanged against any
target; backward-compat aliases preserve the kind16 parity tests.

`pyengine stage` writes `LIB_PREFIX_COMPLEX` to `target_config.cmake`
alongside the existing `LIB_PREFIX`.

#### Results

| Target       | Total tests | mumps tests | Notes                                |
| ------------ | ----------- | ----------- | ------------------------------------ |
| kind16       | 1045        | 23 PASS     | regression unchanged; parity tests now per-target (no behaviour change on kind16) |
| kind10       | (pending)   | 23 build    | builds clean; runtime hang surfaces in fresh stage (separate, pre-existing) |
| multifloats  | (pending)   | 23 build    | builds clean; runtime hits pre-existing migrator bug — MPI_SUM not rewritten to MPI_DD_SUM in mfac_scalings.F (`Allreduce(...USER<contig>, MPI_SUM, ...)` returns "MPI_Op operation not defined for this datatype"). Tracked as B9c |

#### Follow-ups (B9b — RESOLVED 2026-05-02)

- **C-side tests (test_dmumps_c_basic, test_zmumps_c_basic,
  test_dmumps_c_sym) — RESOLVED 2026-05-01 for kind10.** The
  hardcoded `__float128` / `quadmath.h` / `0.0q` literal syntax was
  replaced by a target-aware shim
  (`tests/mumps/c/include/test_real_compat.h`) that picks `__float128`
  for kind16 and `long double` for kind10 along with matching literal
  suffix (`q` / `L`), math functions (`fabsq` / `fabsl`,
  `sqrtq` / `sqrtl`), epsilon (`FLT128_EPSILON` / `LDBL_EPSILON`), and
  snprintf wrapper. The CMake glue (`add_mumps_c_test()`) wires the
  same `TARGET_REAL_HEADER` / `TARGET_REAL_MUMPS_C` /
  `TARGET_REAL_STRUC_C` (plus complex counterparts) macros that
  `c_parity_helpers.c` uses, plus `TEST_TARGET_NAME` and
  `TEST_TARGET_<NAME_UPPER>`.
- **Multifloats C-side tests — RESOLVED 2026-05-02.**
  `test_real_compat.h` gained a `TEST_TARGET_MULTIFLOATS` branch:
  `test_real` is plain `double` (host precision), `test_complex` is
  `{double r, i}`, and `tr_widen` / `tr_narrow` (real) /
  `tc_widen` / `tc_narrow` (complex) translate at the bridge boundary
  to `mumps_float64x2` / `mumps_complex64x2`. The 4×4 hand-built test
  problems don't strain double precision; the C tests verify the
  bridge wires correctly, not multifloats arithmetic in C (Fortran-
  side tests cover that). Drivers got `#ifdef TEST_TARGET_MULTIFLOATS`
  bridge buffers around `id.a` / `id.rhs`.
- **Multifloats Fortran tests — RESOLVED 2026-05-02.** Test drivers
  now `use target_mumps, only: ..., q2t_r, t2q_r` (or `_c` for the
  complex side) and wrap `id%A` / `id%RHS` / `id%A_loc` / `id%RHS_*`
  assignments with `q2t_r(...)` and readbacks with `t2q_r(...)`. The
  `target_conv` library from `tests/common` (already used by lapack /
  blas / scalapack test trees) provides identity converters on
  kind16, intrinsic narrow + widen on kind10, and double-double split
  on multifloats. `target_mumps_body.fypp` re-exports the converters
  through host association so test drivers reach them via their
  existing `use target_mumps` clause. iref_errchk reads of
  `id%RINFOG(6)` (NaN check + bound comparison) cache
  `t2q_r(id%RINFOG(6))` once, replacing the kind16-only
  `real(id%RINFOG(6), kind=ep)` cast.
- **Parity tests target-aware — RESOLVED 2026-05-02.**
  `test_dmumps_c_parity.f90` and `test_zmumps_c_parity.f90` were
  rewritten as `.fypp` and now run on all three targets. The bind(C)
  interface declares the per-target type (`real(16)` / `real(10)` /
  `type(real64x2)`); a per-target staging block converts via
  `q2t_r` / `t2q_r` before/after the C call. CMake processes them
  with `-Dtarget=${TARGET_NAME}`. The kind16-only aliases
  (`c_qmumps_solve` / `c_xmumps_solve`) in `c_parity_helpers.c` are
  dropped — parity tests call the canonical
  `c_real_mumps_solve` / `c_complex_mumps_solve`.
- The `_MUMPS_TESTS_SKIP_EXECUTABLES` skip guard at the head of
  `tests/mumps/CMakeLists.txt` was removed; the kind16-only
  `_c_parity\.f90$` filter was removed; the
  `NOT TARGET_NAME STREQUAL "multifloats"` clause on the C-test loop
  was removed. All three targets now register the same Fortran +
  C-side test sets.

### B9c — Multifloats MUMPS uses `MPI_SUM` on `real64x2` buffers — runtime failure (NEW 2026-05-02)

After Group A unblocked the multifloats mumps test harness (B9b), the
test executables build and link cleanly but fail at runtime with:

```
internal_Allreduce: MPI_Op operation not defined for this datatype
MPI_Allreduce(..., dtype=USER<contig>, MPI_SUM, ...)
```

Root cause: the migrator's `_rewrite_mpi_datatypes` post-pass rewrites
`MPI_DOUBLE_PRECISION → MPI_FLOAT64X2` on the buffer-type argument but
does **not** rewrite the matching `MPI_SUM → MPI_DD_SUM` reduction-op
argument. MPI's built-in `MPI_SUM` is undefined on user-defined
contiguous types (registered via `MPI_Type_contiguous(2, MPI_DOUBLE,
&MPI_FLOAT64X2)` inside `multifloats_mpi_init`), so MUMPS calls like
`MPI_Reduce(..., COUNT, MPI_FLOAT64X2, MPI_SUM, ...)` crash.

Affected files (grep `MPI_SUM` in `/tmp/stg-mf/mumps/src/`):
- `mfac_scalings.F:372,376` — Reduce on real-precision arrays
- `mana_dist_m.F:206,605,626` — Allreduce / scan with real types
- `wsol_distrhs.F`, `msol_distrhs.F` — distributed-RHS reductions
- (others)

Fix locations: `_rewrite_mpi_datatypes` in `src/pyengine/fortran_migrator.py`
should pair the MPI_SUM/MAX/MIN substitution with the buffer-type
substitution. Where `MPI_FLOAT64X2` is the active datatype, use
`MPI_DD_SUM` / `MPI_DD_AMX` / `MPI_DD_AMN`; where `MPI_COMPLEX64X2`
is active, use `MPI_ZZ_SUM` / `MPI_ZZ_AMX` / `MPI_ZZ_AMN`.
`multifloats_mpi.cpp` already exports all six ops via
`multifloats_mpi_init`, and `multifloats_mpi_f.f90` exposes them as
INTEGER handles.

Workaround: none. The test-harness work (B9b) is complete; this bug
blocks runtime success but is independent of the harness rework.

Group A test build verification (2026-05-02):
- kind16: full ctest pass on `/tmp/stage-q` (`23/23` mumps tests)
- multifloats: build clean on `/tmp/stage-mf`; ctest hits B9c on every
  test (datatype + op pairing)
- kind10: build clean on a fresh `/tmp/stg-e-mumps` stage; ctest hangs
  (separate, pre-existing — likely a missing-archive issue in the
  fresh-stage path, not related to Group A changes)
- Pre-existing `_replace_kind_parameter` for multifloats (line 2872
  of fortran_migrator.py) comments out `integer, parameter :: wp =
  kind(1.d0)` lines; `replace_literals` skips literal-wrapping on
  `INTEGER, PARAMETER :: ` lines so degenerate upstream forms like
  `INTEGER, PARAMETER :: ZERO = 0.0D0` (dsol_fwd_aux.F:1093) survive
  the constructor-mode pass intact.
