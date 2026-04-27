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
| UB-06 | pending | `fix/ub-06-reflapack-F-glob` | Add `*.F` to `tests/lapack/reflapack/CMakeLists.txt` glob + EXCLUDE regex | `nm libreflapack_quad.a \| grep iparam2stage_` shows `T` (currently `U`) |
| UB-02 | pending (depends on UB-01) | `fix/ub-02-multifloats-F-param-ordering` | Multifloats `.F` migration places `PARAMETER` lines above `IMPLICIT NONE` | `gfortran -c tsytrd_sb2st.F` clean |
| UB-03 | pending (depends on A+B) | `fix/ub-03-2stage-segfault` OR `doc/ub-03-escalation` | 2-stage segfault in `__GI___libc_free` under `-freal-8-real-16` | `valgrind ctest -R '_2stage'` zero invalid frees |
| UB-04 | pending | `fix/ub-04-cmplx16-locals` | C migrator `cmplx16` alias not applied to local decls; `PB_Cconjg` hardcoded `(double*)` casts | `ctest -R '^pblas_test_pzher2k$'` rel-err ≤ 64·8·k·eps on kind16 |
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
