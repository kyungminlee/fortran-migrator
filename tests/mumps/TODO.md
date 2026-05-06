# tests/mumps — TODO

Resolved items have moved to `CHANGELOG.md`. The original 2026-04-29
"deferred test implementation plan" is closed: tests/CMakeLists.txt
ships, all 26 mumps ctests pass on every target under
`cmake --preset=linux-impi`, and `README.md` documents the run flow.

## D1 — MUMPS 5.8.2 doesn't validate bad inputs (sidestepped, monitoring)

Upstream MUMPS 5.8.2 crashes (SIGSEGV inside analysis or factorization)
on inputs the manual claims are rejected via `INFOG(1) < 0`:

- `id%N = -1` → SIGSEGV before any diagnostic.
- `id%IRN(1) = N + 5` or `id%JCN(1) = N + 3` → SIGSEGV inside
  `qmumps_validate_input_` or callees.
- Mismatched `NNZ` (e.g. `NNZ = 0` with non-empty IRN/JCN/A) → SIGSEGV
  during reformatting.

`target_mumps` exports `check_dmumps_input` / `check_zmumps_input`
plus `MIC_*` codes (see `tests/mumps/common/target_mumps_body.fypp`).
The wrapper layers above MUMPS without touching `id%INFOG(1)`. Tests
at `tests/mumps/fortran/test_{d,z}mumps_errors.f90` exercise BAD_N
(neg + zero), BAD_NNZ, BAD_IRN/JCN (high + zero), SIZE_MISMATCH, plus
a final valid-input pass that reaches MIC_OK and factors via JOB=6.

**Reopen condition**: when MUMPS upstream tightens validation and
returns the documented `INFOG(1) = -16` / `-6` codes, drop the wrapper
and let the tests check `INFOG(1)` directly.

## libmpiseq linkage — kind16/kind10 land, multifloats blocked

Each test source under `tests/mumps/{fortran,c}/` is now built twice
— once linked against real MPI (`mpiexec -n 1`), once linked against
the in-tree `mpiseq` archive (plain binary, suffix `_seq`). libmpiseq
is now Intel-MPI ABI compatible: cmake/CMakeLists.txt overlays the
configured MPI's `mpi.h` / `mpif.h` onto `_mpiseq_src/`, libseq's own
bundled `mpic.c` is replaced by `cmake/mpiseq_c_stubs.c` (compiled
against Intel's signatures), and `pyengine stage` patches libseq's
`mpi.f` to extend `MUMPS_COPY` with `MPI_REAL16` / `MPI_COMPLEX32`
cases the standard upstream dispatch doesn't ship.

**kind16**: 26/26 `_seq` ctests pass; per-test JSON precision reports
are bit-identical to the impi-linked runs (verified via `md5sum`).

**multifloats**: `_seq` binaries STOP at runtime with
`MPI_ALLREDUCE, DATATYPE = 201326592` (Intel's `MPI_DATATYPE_NULL`).
Root cause: the C++ `multifloats_mpi.cpp` registers `MPI_FLOAT64X2`
/ `MPI_DD_SUM` / etc. via `MPI_Type_create_struct` /
`MPI_Op_create`, which our libmpiseq stubs return as `MPI_DATATYPE_NULL`
since they're inert. The runtime handle ends up null, and libseq's
extended `MUMPS_COPY` doesn't recognize null. Fix path: replace the
stubs with a libmpiseq-mode variant that returns specific sentinel
handle values, and extend `MUMPS_COPY`'s dispatch to recognize those
sentinels with the multifloats element sizes (16-byte real, 32-byte
complex). Out of scope for this iteration.

**Reopen condition** (multifloats specifically): an environment that
can't run real MPI needs the multifloats sparse solver, OR a follow-up
adds the multifloats-mode datatype-registration stubs.
