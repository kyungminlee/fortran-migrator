# tests/mumps

Differential precision tests for the migrated MUMPS sparse direct solver.
Mirrors the tests/blas / tests/lapack pattern: every check happens at
`REAL(KIND=16)`, the migrated `${LIB_PREFIX}mumps` archive (`qmumps` for
the kind16 target) is exercised both from Fortran (`call qmumps(id)`) and
from C (`qmumps_c(&id)` via the bridge described in
[`c/include/qmumps_c.h`](c/include/qmumps_c.h)), and per-case JSON
reports are written under `<build>/precision_reports/` for the standard
aggregator.

## How to run

```bash
uv run --project /home/kyungminlee/Code/fortran-migrator \
    python -m pyengine stage /tmp/stg-q --target kind16 \
    --project-root /home/kyungminlee/Code/fm-mumps
cmake -S /tmp/stg-q -B /tmp/stg-q/build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/stg-q/build -j8
ctest --test-dir /tmp/stg-q/build -R '^mumps_' --output-on-failure
```

Tests are wrapped via `mpiexec -n 1` because the migrated qmumps archive
calls MPI primitives unconditionally (see TODO.md B3 for the
libmpiseq-based no-mpiexec alternative — currently a known incomplete
path).

## Layout

```
tests/mumps/
├── CMakeLists.txt        — gates + bridge build + test registration
├── README.md             — this file
├── TODO.md               — open issues, design plan, deferred coverage
├── common/               — shared helpers (prec_kinds, compare,
│                           prec_report, test_data_mumps,
│                           ref_quad_lapack_solve, target_mumps_body.fypp)
├── target_kind16/        — kind16 fypp shim setting prefix=q/x
├── fortran/              — test_*mumps_*.f90 drivers
└── c/
    ├── include/          — quad-precision header overrides for the bridge
    ├── mpiseq_c_stubs.c  — supplementary C-MPI stubs (B3)
    └── test_*mumps_c_*.c — C drivers
```

## Coverage plan

See [TODO.md](TODO.md) for the full design and the per-test parameter
matrix. The first landed tests cover the `JOB=-1 → JOB=6 → JOB=-2`
roundtrip on a small unsymmetric problem; subsequent tests will sweep
SYM, ICNTL ordering, JOB phasing, NRHS, error paths, and the C-side
parity equivalent.
