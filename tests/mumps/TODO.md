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

### UNK-3 — `NO_SAVE_RESTORE` compile flag

The bridge's F77 prototype must agree with the migrated MUMPS's compile-time
`NO_SAVE_RESTORE` setting. The recipe doesn't define the macro, so the
default (SAVE/RESTORE enabled) likely applies — confirm by inspecting the
staged sources for `#if ! defined(NO_SAVE_RESTORE)` blocks and whether
they're active.

### UNK-4 — Link-time satisfaction

Migrated MUMPS declares `depends: [blas, lapack, scalapack]`. Verify that
linking against `${LIB_PREFIX}{mumps,scalapack,lapack,blas}` + `MPI::MPI_Fortran`
is sufficient. The vendored PORD (`external/MUMPS_5.8.2/PORD/`) is C and not
migrated — if the recipe leaves `USE MUMPS_PORD`-callers in the migrated
archive, link will fail with unresolved PORD symbols. Surface during the
build-smoke step.

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

### B7 — Two cmake/recipe trees diverge (fm-mumps vs fortran-migrator)

`/home/kyungminlee/Code/fm-mumps/cmake/CMakeLists.txt` and
`/home/kyungminlee/Code/fortran-migrator/cmake/CMakeLists.txt` are
separate tracked files — same for `recipes/mumps.yaml`.
`pyengine stage` reads from `--project-root` (default:
fortran-migrator). Pass `--project-root /home/kyungminlee/Code/fm-mumps`
to use fm-mumps's copies. Otherwise the two trees drift.

The 2026-04-29 build-pipeline integration ended up touching both:
- `cmake/CMakeLists.txt`: section 9 (MUMPS) added in BOTH trees.
- `recipes/mumps.yaml`: `FAC_FUTURE_NIV2_MOD` added only in fm-mumps.
- `pyengine __main__.py`: `('mumps','mumps.yaml')` added to
  LIBRARY_ORDER (only lives in fortran-migrator).

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

### B4 — Only `kind16` target supported

Recipe `overrides` only lists `kind16` (`mumps_memory_mod_ep.F`,
`mumps_lr_stats_ep.F`). `kind10` / `multifloats` wrappers omitted from
this test directory until upstream support lands.

## Defects discovered during runs

(empty — populated as tests are exercised)
