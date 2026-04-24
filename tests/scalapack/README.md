# Differential precision tests for migrated ScaLAPACK

This directory extends the `tests/pblas/` scaffold to the **distributed**
extended-precision ScaLAPACK libraries (`${LIB_PREFIX}scalapack` —
`qscalapack`, `escalapack`, `ddscalapack`). Each test runs a ScaLAPACK
operation over a BLACS process grid, gathers the distributed result
back to rank 0, and compares against serial quad-precision Netlib
LAPACK applied to the full global matrix.

The structure mirrors `tests/pblas/` closely — the grid/distribution
helpers, per-target wrapper contract, JSON report format, and cyclic
link-group wiring are all inherited. This README only covers what is
*different* from the PBLAS suite.


## Testing model

Identical to PBLAS: every rank materializes the global problem at
REAL(KIND=16) from a shared seed, copies its block-cyclic slab into
a local array, calls `target_pdXXX` (which hands through to the
migrated `pq*` / `px*` entry point), then rank 0 gathers the result
and compares against serial `reflapack_quad` applied to the global
matrix.

Problems chosen to exercise block-cyclic edge cases at small sizes
(`n ∈ {32, 64, 96}`, `mb = nb = 8`) without blowing out runtime.


## Routines covered

| Group          | Routines                                              |
| -------------- | ----------------------------------------------------- |
| linear_solve   | pdgesv, pdgetrf, pdgetrs, pdpotrf, pdpotrs            |
| factorization  | pdgeqrf                                               |
| eigenvalue     | pdsyev (eigenvalues only — `JOBZ='N'`)                |
| svd            | pdgesvd (singular values only)                        |
| auxiliary      | pdlange, pdlacpy, pdlaset                             |

Complex wrappers (`pxgesv`, `pxgeqrf`, `pxheev`) are exposed by the
target module for future expansion; no test drivers are provided yet.


## Reference strategy

Uses `reflapack_quad` (from `tests/lapack/reflapack/`) for the
LAPACK-level comparisons. When `_reflapack_src/` is not staged,
`reflapack_quad` degrades to the system `-llapack` (double precision),
which caps the meaningful precision at ~1e-16 regardless of target —
the tests still run, but the digit count reported will plateau at 15–16.

`refblas_quad` is also linked because the ScaLAPACK↔LAPACK↔BLAS cycle
has to resolve together under `LINK_GROUP(RESCAN)`.


## Dependencies

ScaLAPACK sits on top of the full stack, so every layer below must be
present. The CMake guards skip the suite (with a `STATUS` message) if
any of these is missing:

- `HAVE_REAL16`
- `MPI_Fortran_FOUND`
- `refblas_quad` (from `tests/blas/`)
- `reflapack_quad` (from `tests/lapack/`)
- `${LIB_PREFIX}scalapack`, `${LIB_PREFIX}pblas`, `${LIB_PREFIX}blacs`,
  `${LIB_PREFIX}lapack`
- A per-target wrapper directory `target_${TARGET_NAME}/` with
  `target_scalapack.f90` (currently only `target_kind16/` is populated;
  kind10 and multifloats wrappers are future work)


## Layout

```
tests/scalapack/
├── CMakeLists.txt                    # guards, shim, wrapper, LINK_GROUP cycle
├── common/                           # copied verbatim from tests/pblas/common
│   ├── prec_kinds.f90
│   ├── compare.f90
│   ├── prec_report.f90               # module pblas_prec_report (rank-0 JSON)
│   ├── test_data.f90
│   ├── ref_quad_blas.f90
│   ├── ref_quad_lapack.f90           # copied from tests/lapack/common
│   ├── pblas_grid.f90                # BLACS init + numroc + descinit
│   ├── pblas_distrib.f90             # block-cyclic scatter/gather
│   └── pilaenv.f                     # vendored PBLAS integer helper
├── target_kind16/
│   ├── target_scalapack.f90          # passthrough to pq*/px* routines
│   ├── ptzblas_stubs.f90             # ZZDOTC/ZZDOTU stubs (PTZBLAS gap)
│   └── lamch_shim.f90                # qlamch/qroundup_lwork (INSTALL gap)
├── linear_solve/  test_pdgesv.f90, test_pdgetrf.f90, test_pdgetrs.f90,
│                  test_pdpotrf.f90, test_pdpotrs.f90
├── factorization/ test_pdgeqrf.f90
├── eigenvalue/    test_pdsyev.f90
├── svd/           test_pdgesvd.f90
└── auxiliary/     test_pdlange.f90, test_pdlacpy.f90, test_pdlaset.f90
```

A dedicated `Fortran_MODULE_DIRECTORY`
(`${PROJECT_BINARY_DIR}/fmod/${FORTRAN_MOD_COMPAT_TAG}/scalapack_tests`)
prevents `.mod` collisions with `tests/pblas/` — the helper modules
reuse the same names (`pblas_grid`, `pblas_distrib`, `compare`, …)
because renaming them would defeat the point of reusing the sources.


## Known upstream gaps affecting link

The ScaLAPACK test objects all compile cleanly, but linking the
executables depends on the migrated ScaLAPACK + PBLAS + PBBLAS
archives, which currently have two upstream staging gaps:

1. **`numroc_` / `iceil_` / `ilcm_` missing from PBBLAS closure.**
   These are ScaLAPACK TOOLS integer helpers (precision-independent,
   no renaming needed). `pbxtran*` / `pbxtrnv*` in `libqpbblas` call
   them, but the migrator's recipe does not currently stage the TOOLS
   sources into a linkable archive, so tests whose call graph reaches
   `pbxtran*` fail to link. `test_pdgesv` and `test_pdgetrf` hit this;
   the other 9 drivers happen to avoid those transitive dependencies
   and link cleanly.

2. **`{q,e,dd}lamch` / `roundup_lwork` missing from LAPACK.**
   Same root cause as `tests/lapack`: `external/lapack-3.12.1/INSTALL`
   is declared in the recipe's `extra_symbol_dirs` but not copied by
   `cmd_stage`. Worked around locally by
   `target_${TARGET_NAME}/lamch_shim.f90`, which reconstructs the
   machine-parameter queries from Fortran intrinsics at the target
   precision. Remove the shim once staging covers `INSTALL/`.


## Running

Same as PBLAS — configure a staged target (`uv run python -m pyengine
stage <dir> --target kind16`), build, and run with `mpirun`:

```
cmake -S <stage> -B <stage>/build
cmake --build <stage>/build --target test_pdpotrf -j
ctest --test-dir <stage>/build -R '^scalapack_test_pdpotrf$' --output-on-failure
```

Default `SCALAPACK_TEST_NPROC=4` (2×2 grid). Override at configure time
with `-DSCALAPACK_TEST_NPROC=<n>`.
