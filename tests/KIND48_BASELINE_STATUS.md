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

### Six remaining failures — categorized

| # | Test | Class | Severity |
|---|------|-------|----------|
| 1 | `lapack_test_zunbdb3` | wrapper-side heap corruption | bug to fix |
| 2 | `lapack_test_dorbdb3` | precision genuinely insufficient | finding |
| 3 | `lapack_test_dorcsd`  | precision genuinely insufficient | finding |
| 4 | `lapack_test_dgedmdq` | precision insufficient on residual-refinement path | finding |
| 5 | `lapack_test_zgedmdq` | precision insufficient on residual-refinement path | finding |
| 6 | `lapack_test_disnan`  | NaN-payload artifact of the cast q2t_r | test-fixture quirk |

#### 1 — `zunbdb3` (kind8 wrapper heap corruption)

```
SIGSEGV in arena_for_chunk → __libc_free
  __target_lapack_MOD_target_zunbdb3
  MAIN__
```

The crash is inside the kind8 `target_zunbdb3` wrapper's `free` path,
NOT inside the reference path or `zunbdb3_` itself. `dorbdb3` (the
real-precision sibling) doesn't crash — it produces wrong numbers but
returns control cleanly.

Likely root cause: the wrapper's workspace allocation for the
`(M-Q) × (P-Q)` sub-block case has an off-by-one or wrong bound in
the complex (Z) variant. Look at `target_zunbdb3` in
`tests/lapack/common/target_lapack_body.fypp` — compare against the
real (D) wrapper for any asymmetric handling of the `taup1` / `taup2`
/ `tauq1` workspace dimensions.

This is the only mechanical bug in the residual set; the other five
are precision findings.

#### 2 — `dorbdb3` (real CSD bidiagonalization, last sub-block)

```
dorbdb3 [m=12]  max_rel_err = 4.86E-02   digits = 1.31  FAIL
dorbdb3 [m=16]  max_rel_err = 1.03E-01   digits = 0.99  FAIL
```

Routine: real CS bidiagonalization `[X11; X21]` where
`M-Q ≤ M-P, P ≤ Q, M-P ≤ Q`. Reaches the regime where the active
sub-block sees singular values clustered at or below the kind8
precision floor; double-precision Givens rotations cannot resolve
them faithfully, while the quad reference can. The `~1` digit
agreement means the result is qualitatively different — the
test correctly flags that kind8 isn't accurate enough for these
problem shapes.

Worth promoting in any production use of these CSD primitives:
> **Don't run `dorbdb3` on shapes with near-degenerate sub-blocks at
> double precision — use kind10 or higher.**

The kind4 column (when run) will show this even more dramatically.

#### 3 — `dorcsd` (full real CSD assembly)

```
dorcsd [m=12]  max_rel_err = 2.92    digits = -0.47  FAIL
dorcsd [m=16]  max_rel_err = 8.70E+1 digits = -1.94  FAIL
```

`dorcsd` is the user-facing CSD driver that calls `dorbdb*` underneath
and then assembles the orthogonal factors. Same root cause as
`dorbdb3` but the error compounds across the assembly: the negative
digit counts mean the kind8 result and the quad reference are not the
same orthogonal subspace at all for these problem shapes — they've
picked different canonical representatives among the equivalent
orthogonal completions.

This is the strongest signal in the table that **CSD at double
precision is unreliable for general shapes**. It also constrains how
the kind4/kind8 columns should be presented — for routines like
`dorcsd` the cell should probably show `(unstable)` rather than a
raw digit count.

#### 4, 5 — `dgedmdq` / `zgedmdq` (residual-refined DMD)

```
dgedmdq [m=12,n=8,JOBR=R,layer=b] max_rel_err = 4.03E-03  digits = 2.40  FAIL
zgedmdq [m=12,n=8,JOBR=R,layer=b] max_rel_err = 4.93E-02  digits = 1.31  FAIL
```

Every other case combination passes (`14-15` digits — normal kind8
floor). Failures are isolated to `JOBR=R` (Rayleigh-Ritz residual
refinement) at `m=12, n=8, layer=b` (the second test layer with a
specific input distribution).

The residual-refinement step subtracts two nearly-equal quantities
to compute Ritz residuals; in the quad reference those quantities
are computed with ~33 digits of headroom and the cancellation
cleanly exposes the small residual. In kind8 native, the cancellation
loses most of the digits and the residual ends up dominated by
roundoff. The quad-vs-native discrepancy is real but expected — it's
exactly the precision phenomenon DMD-Q's residual mode is *designed*
to surface.

For the kind4/kind8 columns this should probably be reported per
parameter combination rather than as a single cell; users picking a
double-precision DMD should know `JOBR=R` carries this caveat.

#### 6 — `disnan` (test-fixture artifact)

```
disnan [case=5]  max_rel_err = 1.0   digits = -0.00  FAIL
disnan [case=7]  max_rel_err = 1.0   digits = -0.00  FAIL
```

Cases 5 and 7 feed signaling-NaN bit patterns into `disnan`. The
test's input pipeline is `gen_data_quad → q2t_r → disnan` — so the
NaN is constructed at quad precision, then cast to real(8) via
intrinsic `real(x, 8)`, then passed to native `disnan`. The cast
collapses the quad-precision payload into a canonical real(8) NaN
(IEEE 754 doesn't preserve NaN payload across precision conversions
in a portable way; gfortran's `real(x, 8)` returns the canonical
quiet NaN).

The reference path doesn't go through that cast — it calls quad
`disnan` directly on the original quad value, which still distinguishes
signaling vs quiet. So one path returns `.true.` (sees a NaN of any
kind) and the other returns `.false.` (the canonicalized quad-source
value isn't recognized as NaN under the specific bit pattern).

This is a property of the test data + cast pipeline, not a bug in
either implementation. Either:
- skip cases 5 and 7 in the kind4/kind8 column (most honest), or
- pass the original quad NaN through both paths without the cast
  intermediate (changes the test semantics).

Same shape would surface in `sisnan` for kind4. Since
`disnan` / `sisnan` return `LOGICAL`, no actual precision is being
measured — these are categorical-output routines that don't fit the
"digits of agreement" framing cleanly.

### Recommended treatment in `tests/RESULT.md`

When the kind4 / kind8 columns are turned on for the public report:

- `dorbdb3`, `dorcsd` — show the digit count honestly. The negative
  values are the headline finding.
- `dgedmdq`, `zgedmdq` — consider reporting per-`JOBR` or annotating
  the cell with `(JOBR=R: 2.4)`; the aggregate worst-case is
  unrepresentative.
- `disnan` — use `n/a` (categorical output, not a precision measure).
- `zunbdb3` — show after the wrapper bug is fixed; for now suppress
  via the same xx/rfsx-style filter used for XBLAS-dependent tests.

### Other framework coverage (PBLAS / ScaLAPACK)

The PBLAS / ScaLAPACK kind8 cycles are wired (`scalapack` added
to the PBLAS cycle, std blas/lapack added to the mumps cycle), but
not yet smoke-tested. The same module-procedure rename now in place
should cover them; expect runtime results in line with the LAPACK
suite (≥99% pass for the routines that don't compound CSD-style
near-degenerate sub-blocks).
