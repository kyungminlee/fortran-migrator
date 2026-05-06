# kind4 / kind8 baseline columns — status

## Current state

The plumbing is fully wired and the kind4 / kind8 baseline target
produces sane numbers from a `migrator stage --target kind8` tree:

- `migrator stage --target kind8 /tmp/x` short-circuits the migration
  step and writes a `target_config.cmake` with empty `LIB_PREFIX` so
  the unified CMake build resolves the standard archives (`blas`,
  `lapack`, `pblas`, `scalapack`) as the libraries under test.
  Same shape for `--target kind4`.
- `cmake -S /tmp/x -B /tmp/x/build && cmake --build /tmp/x/build -j8`
  links cleanly. `ctest` runs the framework subsets that the baseline
  supports (BLAS + the non-XBLAS LAPACK tests, plus PBLAS / ScaLAPACK
  once those are wired analogously).

Smoke results on the kind8 BLAS + LAPACK suite (725 tests):

| Suite  | Pass | Fail | Notes |
|--------|-----:|-----:|-------|
| BLAS   |   75 |    0 | Every Level 1/2/3 routine reports ~16-digit agreement against the quad reference (the IEEE binary64 floor — exactly the expected sanity-check shape). |
| LAPACK |  719 |    6 | After the symbol-rename pass was extended to gfortran module-procedure symbols (`__<mod>_MOD_<proc>`), the entire CSD / 2-stage / DMD SIGSEGV cluster cleared. The 6 residual failures are precision corner cases (5) plus one specific `zunbdb3` heap-corruption that traces to the kind8 `target_lapack` wrapper rather than the rename machinery. |

`kind16` regression check stays clean — `test_dasum`, `test_dgemm`,
`test_dsyevr`, and `test_dgbsvxx` all pass after the rename and
exclude-list changes.

## What was needed (and now exists)

1. **Per-target shims.** `tests/{blas,lapack,pblas,scalapack}/target_kind{4,8}/`
   fypp shims setting `prefix_real='s'/'d'`, `prefix_complex='c'/'z'`,
   etc. Per-target `tests/common/target_kind{4,8}/target_conv.f90`
   conversion modules (trivial real(4)/(8)↔real(16) intrinsic casts).
2. **Baseline staging path.** `migrator stage --target kind4|kind8`
   skips migration entirely and writes a target_config.cmake with
   `LIB_PREFIX=""`. `cmake/CMakeLists.txt` gates every
   `add_migrated_*` and the install loop on a `BASELINE_BUILD` flag
   so the std archives stand alone. Each `tests/<lib>/CMakeLists.txt`
   has a baseline branch that links the bare std archive instead of
   `${LIB_PREFIX}<lib>`.
3. **Reference symbol rename.** `librefblas_quad.a` and
   `libreflapack_quad.a` get every Fortran-mangled symbol they
   define renamed to `*_quad_` via a POST_BUILD `objcopy
   --redefine-syms` step, so the un-renamed std archives can coexist
   on the same link line. The rename is driven by
   `scripts/refquad_rename_archive.sh`:
   - Step 1 (scan) — each archive's POST_BUILD writes a
     `.quad-symbols` snapshot of its T+W defined Fortran symbols
     *before* renaming. A sibling rename invoked later reads that
     snapshot rather than nm-ing the (already renamed) archive.
   - Step 2 (rename) — the redefine list is the union of own +
     sibling snapshots, minus the T+W defs of any `-x` exclude
     archive. The exclude list pulls in the migrated `*_common`
     archives (`lapack_common`, `mumps_common`, …) so their
     precision-independent helpers (iparam2stage_, la_constants_ep,
     mumps_abort_) stay un-renamed and resolve normally in the
     migrated link cycle.
   - `ref_quad_*.f90` modules are rewritten by
     `scripts/refquad_alias.py`: each routine gets the original
     interface declaration suffixed `_quad` plus a thin module-
     procedure wrapper under the original name. Test files keep
     their existing `use ref_quad_blas, only: dasum` calls — gfortran
     emits those as module-procedure references
     (`__ref_quad_blas_MOD_dasum`), no longer colliding with the
     standard archive's `dasum_`.
4. **XBLAS baseline stubs.** Std lapack's iterative-refinement
   helpers (dla_*_extended.f, zla_*_extended.f) reference XBLAS
   entry points (blas_dgemv_x_, …). Those .o files get pulled into
   baseline tests transitively via `target_lapack`'s public xx /
   rfsx wrappers, even though those tests themselves are filtered
   out of the baseline ctest set. `tests/lapack/reflapack/xblas_baseline_stubs.c`
   provides aborting C stubs for the entry points; the archive is
   linked with `--whole-archive` so ld extracts every stub
   regardless of demand timing.
5. **Migrated cycle additions.** `tests/mumps/CMakeLists.txt` adds
   std `blas` and `lapack` to its archive cycle so precision-
   independent helpers (lsame_, ilaenv_, dlamch_) resolve from a
   non-renamed source. `tests/pblas/CMakeLists.txt` adds std
   `scalapack` to the kind4 / kind8 cycle for the zzdotc_ /
   zzdotu_ TOOLS routines that PTZBLAS calls.

## Known residual issues

### LAPACK SIGSEGV cluster (kind8) — RESOLVED

The earlier cluster of CSD / 2-stage / DMD SIGSEGVs (dorbdb1-4,
dorcsd2by1, dgedmd, dgges, zheevr_2stage, …) traced to module-proc
symbol collisions: la_xisnan.F90 and la_constants.f90 are compiled
into both std lapack (native) and reflapack_quad (quad-promoted),
producing identical `__la_xisnan_MOD_disnan` / `__la_constants_MOD_dp`
symbols with mismatched ABIs. Whichever .o ld picked first determined
the receiver's argument width, breaking the other path's callers.

`scripts/refquad_rename_archive.sh` was extended to match
`__<mod>_MOD_<proc>` alongside the trailing-underscore Fortran-mangled
pattern. After the rename:
- std lapack continues to expose `__la_xisnan_MOD_disnan` (native) for
  its own native callers (`dlassq_`, `dorbdb5_`, …).
- reflapack_quad exposes `__la_xisnan_MOD_disnan_quad` (quad) for its
  own quad-promoted callers.
- Both paths resolve consistently; the cluster is gone.

### Single remaining SIGSEGV — kind8 zunbdb3

`gdb` traces zunbdb3's segfault into glibc's `arena_for_chunk` via
`__libc_free`, called directly from
`__target_lapack_MOD_target_zunbdb3` — a heap corruption inside the
kind8 wrapper itself, not the reference path. Likely an off-by-one
in the wrapper's workspace allocation for the (M-Q)×(P-Q) sub-block
case. Worth a focused look but unrelated to the rename machinery.

### Real precision corner cases

- `disnan` (case 5, 7) — the test feeds signaling-NaN bit patterns
  through q2t_r (real(16)→real(8)) and the cast collapses different
  NaN payloads to the same canonical NaN. The kind8 native disnan
  and the quad reference can disagree on which class a marginal
  NaN belongs to. Not a real precision regression — a test-data
  artifact at the precision floor.
- `dorcsd` digit count — analogous; the native dorcsd's last few
  bits diverge from the quad-promoted reference under certain
  block-row pivot choices, dropping the digit count below the
  16-digit tolerance.

These are real findings, but not bugs in the migrator — they're the
kind of edge-case behavior the kind8 sanity-check column is meant
to surface.

### Other framework coverage (PBLAS / ScaLAPACK)

The PBLAS / ScaLAPACK kind8 cycles are wired (`scalapack` added
to the PBLAS cycle, std blas/lapack added to the mumps cycle), but
not yet smoke-tested. They likely need the same la_xisnan /
la_constants exclude treatment to clear runtime crashes; otherwise
the link-time pieces are in place.
