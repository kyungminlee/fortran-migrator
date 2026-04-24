# LAPACK differential precision tests

Differential precision tests for the migrated LAPACK libraries
(`qlapack` / `elapack` / `ddlapack`), following the scheme described in
`tests/README.md`. Each test runs a routine through the migrated
extended-precision implementation and through a REAL(KIND=16) reference
(`reflapack_quad`), then compares the outputs in KIND=16 and writes a
per-test JSON report.

Routines covered (15):

| group          | routines                                                |
|----------------|---------------------------------------------------------|
| linear_solve   | dgesv, dgetrf, dgetrs, dpotrf, dpotrs, zgesv            |
| factorization  | dgeqrf, dorgqr, zgeqrf                                  |
| eigenvalue     | dsyev, zheev                                            |
| svd            | dgesvd                                                  |
| auxiliary      | dlange, dlacpy, dlaset                                  |

For `dsyev`, `zheev`, and `dgesvd` the eigenvector / singular-vector
columns have a precision-sensitive sign convention (an internal LAPACK
test on a near-zero value can flip across precisions). The tests
therefore compare only the canonical unique outputs (eigenvalues W,
singular values S) directly, and validate the vectors through the
sign-immune invariants `||A·Z − Z·diag(W)||/||A||` and
`||A − U·diag(S)·Vᵀ||/||A||`.


## Known migrator gaps

Two pre-existing gaps in the migrator pipeline affect these tests. Both
are worked around locally within `tests/lapack/`; the workarounds
become redundant and can be removed once the pipeline fixes land.

### Gap 1 — `cmd_stage` does not stage LAPACK reference sources

**Symptom:** `pyengine stage <dir> --libraries lapack` does not copy
`external/lapack-3.12.1/SRC/` into the staging tree. `tests/blas` has
an analogous copy of `BLAS/SRC/` to `_refblas_src/`, but no such copy
is made for LAPACK.

**Effect:** `reflapack/CMakeLists.txt` looks for `_reflapack_src/` at
the staging root and, when absent, falls back to linking the system
`-llapack` (double precision). For `kind16` / `multifloats` targets the
explicit REAL(KIND=16) interfaces in `common/ref_quad_lapack.f90` then
pass 16-byte arguments to a KIND=8 reference routine — an ABI
mismatch that produces garbage reference values. Every test will
fail. For a `kind10` target the mismatch is the same but the damage
is less visible.

**Fix location:** `src/pyengine/__main__.py`, `cmd_stage`. The existing
BLAS copy is around the staging of `_refblas_src/`; add an analogous
copy of `external/lapack-3.12.1/SRC/*.{f,f90,F90}` to
`_reflapack_src/`, plus `external/lapack-3.12.1/INSTALL/dlamch.f` and
`droundup_lwork.f` (used by many SRC routines).

**Local workaround for verification:** from the repo root,

```bash
python -m pyengine stage /tmp/staging --target kind16 --libraries blas lapack
mkdir -p /tmp/staging/_reflapack_src
cp external/lapack-3.12.1/SRC/*.f   /tmp/staging/_reflapack_src/
cp external/lapack-3.12.1/SRC/*.f90 /tmp/staging/_reflapack_src/
cp external/lapack-3.12.1/SRC/*.F90 /tmp/staging/_reflapack_src/
cp external/lapack-3.12.1/INSTALL/dlamch.f           /tmp/staging/_reflapack_src/
cp external/lapack-3.12.1/INSTALL/droundup_lwork.f   /tmp/staging/_reflapack_src/
cmake -S /tmp/staging -B /tmp/staging/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/staging/build -j8
ctest --test-dir /tmp/staging/build -R '^lapack_'
```

Repeat with `--target kind10` and `--target multifloats` to verify all
three targets.


### Gap 2 — DLAMCH / DROUNDUP_LWORK are never migrated

**Symptom:** `recipes/lapack.yaml` declares
`extra_symbol_dirs: [external/lapack-3.12.1/INSTALL]`. That directive
makes the symbol scanner aware of `dlamch`, `droundup_lwork`, etc. for
the rename map, but the migration loop in `cmd_stage` does not
actually copy and rewrite those files. As a result `qlamch`, `elamch`,
`ddlamch`, and the matching `*roundup_lwork` symbols are missing from
the migrated library.

**Effect:** The migrated lapack library links with unresolved
references to `qlamch_` (called by `qlascl`, `qsyev`, `qgesvd`,
`qgetrf2`, …), and every LAPACK test fails to link.

**Fix location:** extend the migration loop in `src/pyengine/__main__.py`
(or the relevant pipeline stage) to include files discovered via
`extra_symbol_dirs`, not just the recipe's `source_dir`.

**Local workaround:** each `target_<target>/` directory ships a
`lamch_shim.f90` that reconstructs `{q,e,dd}lamch` and
`{q,e,dd}roundup_lwork` using Fortran intrinsics (for kind10/kind16)
or the well-known double-double identities (for multifloats). The
shim is compiled into a separate static library,
`lapack_test_shim`, which `tests/lapack/CMakeLists.txt` links into
the test executables *after* `${LIB_PREFIX}lapack`. Single-pass
linkers therefore pull the shim only for symbols the main library
leaves unresolved — and once the migrator is fixed and the real
`qlamch` appears in `${LIB_PREFIX}lapack`, the shim silently stops
being pulled in. Delete the `lamch_shim.f90` files, along with the
`lapack_test_shim` block in `tests/lapack/CMakeLists.txt`, once the
migrator gap is closed.


## Retiring the workarounds

After both gaps are fixed, the following changes to `tests/lapack/`
collapse it back to the BLAS scaffold shape:

1. Delete `target_kind10/lamch_shim.f90`, `target_kind16/lamch_shim.f90`,
   `target_multifloats/lamch_shim.f90`.
2. Delete the `# ── LAMCH shim …` block (the `lapack_test_shim`
   library) at the bottom of `tests/lapack/CMakeLists.txt`.
3. Drop the "NOTE: `pyengine stage` does not yet copy LAPACK sources…"
   paragraph from `tests/lapack/reflapack/CMakeLists.txt`.

Everything else (the wrappers, common modules, tests) stays
unchanged.
