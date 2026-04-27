# Test TODO drain — upstream migrator items

This doc tracks the upstream-migrator project that drains the bugs
documented in the various `tests/<lib>/TODO.md` files. The first cut is
**upstream-only**: each item is a defect in `src/pyengine/`, a recipe, or
build infra that surfaced while the differential precision tests were
being written. Library-local test additions (the `LG-*` items in the
backing plan) are out of scope.

Backing plan: `~/.claude/plans/start-a-project-to-stateless-bumblebee.md`.

## Status legend

- `pending` — not started
- `in-progress` — branch open, fix being implemented
- `merged` — fix landed on `develop`
- `escalated` — root cause is upstream of fortran-migrator; tracked in `doc/escalations/`

## Items

| ID | Status | Branch | One-line | E2E verify |
|----|--------|--------|----------|-----------|
| UB-01 | merged | `fix/ub-01-F-glob` | Add `*.F` to PyEngine source globs in `cmd_stage`/`cmd_verify` | `nm liblapack_common-*.a \| grep iparam2stage_` shows `T` (defined) |
| UB-06 | merged | `fix/ub-06-reflapack-F-glob` | Add `*.F` to `tests/lapack/reflapack/CMakeLists.txt` glob + EXCLUDE regex | `nm libreflapack_quad.a \| grep iparam2stage_` shows `T` (defined) |
| UB-02 | merged | `fix/ub-02-multifloats-F-param-ordering` | Multifloats `.F` migration places `PARAMETER` lines above `IMPLICIT NONE` | `gfortran -c tsytrd_sb2st.F` clean |
| UB-03 | pending (depends on A+B) | `fix/ub-03-2stage-segfault` OR `doc/ub-03-escalation` | 2-stage segfault in `__GI___libc_free` under `-freal-8-real-16` | `valgrind ctest -R '_2stage'` zero invalid frees |
| UB-04 | merged | `fix/ub-04-cmplx16-locals` | C migrator `cmplx16` alias not applied to local decls; `PB_Cconjg` hardcoded `(double*)` casts | `ctest -R '^pblas_test_pzher2k$'` rel-err ≤ 64·8·k·eps on kind16 |
| UB-05 | pending | `fix/ub-05-extern-c-body-wrap` | Multifloats PBLAS `extern "C"` doesn't propagate to per-file `<name>_.c` defs | `nm libqpblas.a \| grep ' T pdgemm_'` un-mangled |

## UB-01 — merged

**What**: `Path.glob` is case-sensitive on Linux, so `src_dir.glob('*.f')`
in `cmd_stage`/`cmd_verify` did not match `iparam2stage.F`,
`dsytrd_sb2st.F`, or `zhetrd_hb2st.F`. The migrator wrote the files into
the staging dir but the manifest never picked them up, so linking
against the migrated LAPACK failed with `undefined reference to
iparam2stage_`.

**How**: extracted `_collect_source_files(src_dir, language)` helper
covering all four cases (`*.f` / `*.F` / `*.f90` / `*.F90`), with
inode-based dedupe so case-insensitive filesystems do not double-stage.
`cmd_stage` and `cmd_verify` now call it. Unit tests in
`src/tests/test_main_globs.py`.

**Verify**:

```bash
uv run pytest src/tests/test_main_globs.py -v
rm -rf /tmp/stg-k16-ub01
uv run python -m pyengine stage /tmp/stg-k16-ub01 --target kind16 --libraries blas lapack
cmake -S /tmp/stg-k16-ub01 -B /tmp/stg-k16-ub01/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg-k16-ub01/build -j8
nm /tmp/stg-k16-ub01/build/liblapack_common-*.a | grep iparam2stage_  # T iparam2stage_
nm /tmp/stg-k16-ub01/build/libqlapack-*.a | grep qsytrd_sb2st_       # T qsytrd_sb2st_
```

**Note**: with UB-01 merged, the next blocker for any 2-stage runtime
test is UB-06 — `libreflapack_quad.a` still has undefined references
because `tests/lapack/reflapack/CMakeLists.txt` has the same glob bug.

## UB-06 — merged

**What**: `tests/lapack/reflapack/CMakeLists.txt` enumerated reference
sources via `file(GLOB ${_reflap_dir}/*.f)` / `*.f90` / `*.F90` but not
`*.F`, so the quad-promoted Netlib reference library
`libreflapack_quad.a` lacked `iparam2stage`, `dsytrd_sb2st`,
`ssytrd_sb2st`, `chetrd_hb2st`, `zhetrd_hb2st`. Even with UB-01 fixed,
linking any 2-stage test against the reference library failed with
`undefined reference`.

**How**: added `file(GLOB _reflap_F ${_reflap_dir}/*.F)` and renamed
`_reflap_FF` → `_reflap_F90` for clarity. Widened the EXCLUDE REGEXes
that filter mixed-precision and helper files from `\.f9?0?$` to
`\.[fF]9?0?$` so a hypothetical `*.F` rename of a bridge file would
still be filtered correctly.

**Verify**:

```bash
rm -rf /tmp/stg-k16-ub06
uv run python -m pyengine stage /tmp/stg-k16-ub06 --target kind16 --libraries blas lapack
cmake -S /tmp/stg-k16-ub06 -B /tmp/stg-k16-ub06/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg-k16-ub06/build -j8 --target reflapack_quad
nm /tmp/stg-k16-ub06/build/tests/lapack/reflapack/libreflapack_quad.a \
  | awk '/^[0-9a-f]+ T (iparam2stage_|dsytrd_sb2st_|ssytrd_sb2st_|chetrd_hb2st_|zhetrd_hb2st_)$/'
# All 5 symbols print as T (defined).
```

**Note**: with UB-01 + UB-06 merged, 2-stage tests can now LINK. The
remaining blockers are UB-02 (multifloats `.F` PARAMETER ordering for
the migrated side) and UB-03 (quad-promoted reference segfault at
runtime). Either can be tackled next.

## UB-02 — merged

**What**: when migrating capital-F preprocessed Fortran (e.g.
`dsytrd_sb2st.F`, `zhetrd_hb2st.F`) to the multifloats target, the
PARAMETER-conversion anchor walker in `insert_use_multifloats` stopped
at the `#if defined(_OPENMP)` directive sitting between the USE clause
and the IMPLICIT NONE line. It treated `#if` as the first executable
statement and inserted the converted `PARAMETER` assignments above
`IMPLICIT NONE`. gfortran then rejected the file with `Unexpected
IMPLICIT NONE statement`.

**How**: extended the decl-block walker in
`src/pyengine/fortran_migrator.py` (the `insert_use_multifloats`
function) to skip lines starting with `#`. The body of the `#if/#endif`
block in this position is itself decl-only (`use omp_lib`) and is
already matched by the existing `USE` keyword check, so just skipping
the directive markers is sufficient. Regression test
`test_param_assignment_lands_below_implicit_none_with_openmp_directive`
in `src/tests/test_multifloats_transforms.py` covers a synthetic
fixture with the same OPENMP-guard pattern and asserts ordering.

**Verify**:

```bash
uv run pytest src/tests/test_multifloats_transforms.py -k 'openmp_directive' -v
rm -rf /tmp/stg-mf-ub02
uv run python -m pyengine stage /tmp/stg-mf-ub02 --target multifloats --libraries blas lapack
cmake -S /tmp/stg-mf-ub02 -B /tmp/stg-mf-ub02/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg-mf-ub02/build -j8 --target tlapack
nm /tmp/stg-mf-ub02/build/libtlapack-*.a | awk '/^[0-9a-f]+ T (tsytrd_sb2st_|vhetrd_hb2st_)$/'
nm /tmp/stg-mf-ub02/build/liblapack_common-*.a | awk '/^[0-9a-f]+ T iparam2stage_$/'
# All three symbols print as T (defined).
```

**Note**: with UB-01 + UB-02 + UB-06 merged, 2-stage routines now
compile and link on all three targets. The remaining blocker for
runtime is UB-03 (quad-promoted reference segfault under
`-freal-8-real-16`).

## UB-04 — merged

**What**: PBLAS `c_type_aliases` (`cmplx16` → `cmplxQ`/`cmplxE`/`cmplxDD`)
was applied only when c_migrator CLONED a precision-family member
(via `_apply_c_type_subs`). Precision-independent dispatchers
(`PB_Ctzhemm.c`, `PB_Ctzher2k.c`, `PB_Cconjg.c`, …) are flat-copied
through the copy-originals path which did not run any alias pass —
so they kept their `cmplx16 Calph16;` local declarations and their
hardcoded `((double*)CALPHA)[REAL_PART]` strides. On kind16 the
buffer (`cmplxQ`) is 32 bytes but `(double*)`-strided writes only
copy 8 of those bytes, leaving the high half garbage. The
`pzher2k` rank-2k test surfaced this as `max-rel-err ≈ 1.0`.

**How**: added a new helper `_apply_aliases_to_original` in
`src/pyengine/c_migrator.py` that applies the recipe's
`c_type_aliases` rules AND a new `c_pointer_cast_aliases` rule (also
new, plumbed through `src/pyengine/config.py` and the pipeline) to
copy-original C sources. Distinct from the cloned-file pass because
it omits the broad `double` → `REAL_TYPE` substitution — copy-originals
contain precision-dispatch logic (`switch( TYPE->type ) { case DCPLX:
… }`) where the bare `double` keyword identifies the dispatch arm
and must stay. Only the casts and aliased typedefs need rewriting.

`recipes/pblas.yaml` declares the new cast rule for `(double*)` and
`(float*)` → `({REAL_TYPE}*)`. The bug also folds in **LG-03**
(re-verify pzher2k after fix) — verification of UB-04 IS LG-03's
verification.

**Verify**:

```bash
uv run pytest src/tests/test_c_migrator_multifloats.py -k 'aliases_to_original' -v
rm -rf /tmp/stg-k16-ub04
uv run python -m pyengine stage /tmp/stg-k16-ub04 --target kind16 \
  --libraries blas blacs lapack ptzblas pbblas pblas
grep cmplx16 /tmp/stg-k16-ub04/pblas/src/PB_Ctzher2k.c             # zero hits
grep '(double\*)' /tmp/stg-k16-ub04/pblas/src/PB_Cconjg.c          # zero hits
cmake -S /tmp/stg-k16-ub04 -B /tmp/stg-k16-ub04/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg-k16-ub04/build -j8
ctest --test-dir /tmp/stg-k16-ub04/build --output-on-failure
# 109/109 ctest tests pass on develop branch (pzher2k test ships only
# on the `tests` branch and goes from red→green there once `tests`
# rebases onto develop).
```

**Note**: this fix touches recipe-driven C migration only — multifloats
PBLAS still requires UB-05 (`extern "C"` body wrap) before its
Fortran-callable entry points become reachable.
