# Differential precision tests for migrated PBLAS

This directory extends the `tests/blas/` scaffold to the **distributed**
extended-precision PBLAS libraries (`${LIB_PREFIX}pblas` —
`qpblas`, `epblas`, `ddpblas`). Each test runs a PBLAS operation over
a BLACS process grid, gathers the distributed result back to rank 0,
and compares against serial quad-precision BLAS applied to the full
global matrix.

The structure mirrors `tests/blas/` closely — the core testing model
(reference vs. target, tolerance formula, JSON report, per-target
wrapper contract) is unchanged. Everything below describes only the
things that are *different* because the routines under test are
distributed.


## Testing model

Each test runs in the following order (on every MPI rank):

1. `grid_init()` — BLACS + MPI initialization, picks a 2D grid shape
   as close to square as possible for the given process count
   (2×2 for NP=4, 1×1 for NP=1, 1×N for prime NP).
2. For each problem shape:
   - Every rank independently materializes the **full global**
     vector/matrix at REAL(KIND=16) using a deterministic seed — so
     ranks don't need to scatter, they each read their own
     block-cyclic piece from their own mirrored copy.
   - The per-target wrapper (`target_pdgemm` etc.) converts the local
     arrays to the target's precision, calls the migrated PBLAS
     routine (`pqgemm` / `pegemm` / `pddgemm`), and converts the
     in/out buffers back to quad.
   - Rank 0 gathers the distributed result from every process via
     point-to-point `MPI_BYTE` transfers (MPI_REAL16 is not portable).
   - Rank 0 computes the reference by running the serial
     **quad-precision Netlib BLAS** (`refblas_quad`, built by
     `tests/blas/refblas/`) on the full global matrix.
   - Rank 0 computes `max_rel_err` and calls `report_case`; non-zero
     ranks skip it.
3. `report_finalize()` is collective — rank 0 broadcasts the pass/fail
   flag so every rank exits with the same code and ctest observes a
   uniform failure/success across the entire `mpirun` invocation.
4. `grid_exit()` tears down BLACS + MPI.

Tolerance still follows `safety_factor * n_ops * target_eps`, but the
safety factor is doubled over the BLAS suite (16 → 32 for most cases,
64 for triangular solves) to absorb the small amount of extra rounding
that comes from SUMMA-style summation orders in PBLAS kernels.


### Reference strategy

PBLAS is a C library, so the `-freal-8-real-16` promotion that builds
`refblas_quad` from Netlib Fortran sources does not apply. Rather than
porting 122 C files to `__float128` plus MPI types, this test suite
uses the **serial** quad-precision BLAS as a differential reference:
the migrated distributed routine at target precision should agree with
serial BLAS at quad precision up to rounding and SUMMA summation order.

Consequence: `tests/blas/` must also be staged for PBLAS tests to
build — the CMakeLists gates on `TARGET refblas_quad` and skips
gracefully otherwise.


## Directory layout

```
tests/pblas/
├── README.md                     ← you are here
├── CMakeLists.txt                ← gated build, registers ctest entries
├── common/
│   ├── prec_kinds.f90            ← ep=16, dp=8 (shared with tests/blas)
│   ├── compare.f90               ← rel-error helpers (shared)
│   ├── test_data.f90             ← seeded RNG (shared — same seed → same
│   │                               data on every rank, no scatter needed)
│   ├── prec_report.f90           ← rank-0-only JSON writer
│   ├── ref_quad_blas.f90         ← subset of serial quad-BLAS interfaces
│   ├── pblas_grid.f90            ← BLACS init/exit + grid shape picker
│   │                               + inline numroc / descinit / g2l
│   └── pblas_distrib.f90         ← block-cyclic scatter (mirror → local)
│                                   + gather (local → rank-0 global)
├── target_kind10/target_pblas.f90
├── target_kind16/target_pblas.f90
├── target_multifloats/target_pblas.f90
├── level1/                       ← 8 Level 1 PBLAS tests
│   ├── test_pddot.f90            ← scalar out
│   ├── test_pdnrm2.f90
│   ├── test_pdasum.f90
│   ├── test_pdscal.f90           ← in-place vector
│   ├── test_pdaxpy.f90
│   ├── test_pdcopy.f90           ← out-of-place vector
│   ├── test_pzdotc.f90           ← complex
│   └── test_pzaxpy.f90
├── level2/                       ← 5 Level 2 PBLAS tests
│   ├── test_pdgemv.f90
│   ├── test_pdger.f90            ← rank-1 update
│   ├── test_pdsymv.f90
│   ├── test_pdtrsv.f90           ← triangular solve (diag-dominant)
│   └── test_pzgemv.f90
└── level3/                       ← 7 Level 3 PBLAS tests
    ├── test_pdgemm.f90
    ├── test_pdsymm.f90
    ├── test_pdsyrk.f90
    ├── test_pdtrmm.f90
    ├── test_pdtrsm.f90
    ├── test_pzgemm.f90
    └── test_pzherk.f90
```


## Per-target wrapper contract

Every `target_<target>/target_pblas.f90` exposes the same module
`target_pblas` with:

```fortran
module target_pblas
    public :: target_name, target_eps
    public :: target_pddot, target_pdgemm, …  ! one per routine
end module
```

- `target_name`: used in JSON filename (`<routine>.<target_name>.json`).
- `target_eps`: target's machine epsilon for tolerance sizing.
- Each wrapper takes REAL(KIND=16) / COMPLEX(KIND=16) local
  distributed arrays + unchanged descriptor/index arguments, casts
  element-wise to the target native type, calls the migrated PBLAS
  routine, and casts the in/out buffer back to quad. Descriptors are
  distribution info — they pass through unchanged regardless of target
  precision.

The migrated name scheme (scalapack-style single-character slot
substitution):

| original | kind10   | kind16   | multifloats |
|----------|----------|----------|-------------|
| `pdgemm` | `pegemm` | `pqgemm` | `pddgemm`   |
| `pddot`  | `pedot`  | `pqdot`  | `pdddot`    |
| `pzgemm` | `pygemm` | `pxgemm` | `pzzgemm`   |
| `pzherk` | `pyherk` | `pxherk` | `pzzherk`   |


## Running

```bash
python -m pyengine stage /tmp/staging-kind16 --target kind16 \
       --libraries blas blacs ptzblas pbblas pblas
cmake -S /tmp/staging-kind16 -B /tmp/staging-kind16/build \
      -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/staging-kind16/build -j8
ctest --test-dir /tmp/staging-kind16/build -R '^pblas_' --output-on-failure
python scripts/precision_report.py /tmp/staging-kind16/build/precision_reports/
```

The process count is fixed at configure time via `-DPBLAS_TEST_NPROC=<N>`
(default 4 → 2×2 grid). `-DPBLAS_TEST_NPROC=1` runs every test in a
degenerate 1×1 grid, useful as a sanity check.

Staging **must** include `blas`, `blacs`, `ptzblas`, `pbblas`, and
`pblas`. `cmake/CMakeLists.txt` wires `${LIB_PREFIX}pblas` to depend
only on `blas` + `blacs`, but the migrated pblas C code calls into
PTZBLAS (e.g. `qvasum_`) and PBBLAS kernels at runtime — the tests
link those two libraries explicitly to satisfy the references. Without
the `blas` library staged, `refblas_quad` isn't built and the PBLAS
tests skip.


## Reading the report

The JSON output uses the same filename convention as the BLAS suite
(`<routine>.<target>.json` in `<build>/precision_reports/`) and
`scripts/precision_report.py` scans that directory without any
PBLAS-specific changes — PBLAS rows appear in the same summary table
as BLAS rows.

| routine | kind10 | kind16 | multifloats |
|---------|--------|--------|-------------|
| ddot    | 17.59  | exact  | 30.68       |
| pddot   | ~17.2  | exact  | ~29.8       |
| pdgemm  | ~17.0  | exact  | ~29.5       |

For kind16 PBLAS rows, "exact" is the expected outcome when the
tolerance comfortably absorbs SUMMA summation-order differences —
reference and target are both at REAL(KIND=16), and the migrated
distributed algorithm differs from serial BLAS only in associativity.


## Adding a new test

For a PBLAS routine `pxroutine` already wrapped in every
`target_<target>/target_pblas.f90` and declared in
`common/ref_quad_blas.f90`:

1. Drop a new program into the appropriate `level{1,2,3}/` directory
   (one of the existing `test_pd*.f90` files makes a suitable
   starting template).
2. `CMakeLists.txt` picks up `test_*.f90` via
   `file(GLOB CONFIGURE_DEPENDS …)`, so no manual registration is
   needed.

For a routine **not** yet declared:
1. Add an explicit interface to `common/ref_quad_blas.f90` for the
   corresponding serial reference routine.
2. Add wrappers to *every* `target_<target>/target_pblas.f90` —
   signature is uniformly REAL(KIND=16) / COMPLEX(KIND=16) on the
   outside.


## Upstream issues to escalate

Bringing up the first test that actually links `${LIB_PREFIX}pblas`
into an executable surfaced five gaps in the migrator / build pipeline.
Every gap is papered over locally under `tests/pblas/` so that the
test suite builds; each workaround should be replaced by a proper fix
upstream, at which point the corresponding hack here can be deleted.

### 1. BLACS doesn't expose its C API (`Cblacs_*`)

**Symptom.** Linking any executable against `${LIB_PREFIX}pblas`
fails with undefined references to `Cblacs_gridinfo`, `Cblacs_pinfo`,
`Cblacs_get`, `Cblacs_gridexit`, `Cblacs_abort`, `Cblacs_pnum`,
`Csys2blacs_handle`, `C{q,e,dd}gsum2d`, `C{q,e,dd}gamx2d`, etc.

**Root cause.** Upstream BLACS sources in
`external/scalapack-2.2.3/BLACS/SRC/` use `#if (INTFACE == C_CALL)` to
gate between the Fortran-callable (`blacs_gridinfo_`) and the
C-callable (`Cblacs_gridinfo`) entry points. Upstream's
`BLACS/SRC/CMakeLists.txt` compiles every support + communication
source *twice* — once unmodified, once with `-DCallFromC` via a
`src-C.c.in` template — and combines the two object sets into a
single `libblacs.a`. The migrator's `add_migrated_c_library(blacs)`
compiles the sources only once, without `-DCallFromC`, so the
resulting `${LIB_PREFIX}blacs` / `blacs_common` archives only expose
the Fortran-callable half. PBLAS C code is written against the C
half, so every downstream link of PBLAS breaks.

**Local workaround.** `tests/pblas/CMakeLists.txt` re-compiles the
dual-interface BLACS sources (`blacs_*.c` + `{sys2blacs,blacs2sys,
free_handle}_.c` + the 11 `*{g,t}*2d_.c` communication families) with
`-DCallFromC` into a `pblas_test_blacs_c_api` static library, and adds
it to the link cycle.

**Fix upstream.** Replicate the upstream dual-compile trick inside
`add_migrated_c_library` (or specifically for BLACS) so the generated
BLACS library ships both interfaces. Once that lands, delete the
`pblas_test_blacs_c_api` target from `tests/pblas/CMakeLists.txt`.

### 2. `${LIB_PREFIX}pblas` link graph missing PTZBLAS / PBBLAS

**Symptom.** Linking yields `undefined reference to qvasum_`,
`pbqmatadd_`, and similar serial-BLAS-like kernels.

**Root cause.** `cmake/CMakeLists.txt` wires `${LIB_PREFIX}pblas`'s
`target_link_libraries` to `blas` + `blacs` only, but the migrated
PBLAS C entry points call into serial kernels that live in
`${LIB_PREFIX}ptzblas` (e.g. `qvasum_`, `qvvdotc_`) and
`${LIB_PREFIX}pbblas` (e.g. `pbqmatadd_`, `pbqtran_`). Users are
currently expected to know to link those libraries themselves.

**Local workaround.** `pblas_test_target` in
`tests/pblas/CMakeLists.txt` adds `${LIB_PREFIX}ptzblas` and
`${LIB_PREFIX}pbblas` to its PUBLIC link set when those targets exist.

**Fix upstream.** Extend the dependency block around line 299 of
`cmake/CMakeLists.txt`:
```cmake
if(TARGET ${LIB_PREFIX}pblas)
    foreach(_dep ${LIB_PREFIX}blas ${LIB_PREFIX}blacs
                 ${LIB_PREFIX}ptzblas ${LIB_PREFIX}pbblas)
        ...
```

### 3. `pilaenv.f` is unreachable via the PBLAS recipe

**Symptom.** `undefined reference to pilaenv_` from PBLAS C helpers
(`PB_Cplascal`, `PB_Cpsymm*`, `PB_Cptrmm*`, `PB_Cptrsm*`, etc.).

**Root cause.** `recipes/pblas.yaml` sets `extensions: [.c, .h]`,
which filters out the sole Fortran file in `PBLAS/SRC/` — `pilaenv.f`,
a precision-independent block-size query. The file is still referenced
by PBLAS C code at link time but never compiled into any migrated
library.

**Local workaround.** A verbatim copy of `pilaenv.f` lives at
`tests/pblas/common/pilaenv.f` and is compiled into a
`pblas_test_pilaenv` static library.

**Fix upstream.** Either (a) widen `recipes/pblas.yaml`'s `extensions`
to include `.f`, and teach the recipe to copy `pilaenv.f` into the
staged `pblas/src/` tree; or (b) call out `pilaenv.f` as a non-renamed
precision-independent source that gets directly added to the
`${LIB_PREFIX}pblas` build via `add_migrated_c_library`.

### 4. Dangling `ZZDOTC` / `ZZDOTU` in migrated PTZBLAS

**Symptom.** Linking kind10 or kind16 PBLAS fails with
`undefined reference to zzdotc_` / `zzdotu_`, pulled in from
`[xy]vvdotc_.f` in `${LIB_PREFIX}ptzblas`.

**Root cause.** Upstream `PBLAS/SRC/PTZBLAS/zvvdotc.f` calls
`ZZDOTC` — a ScaLAPACK-TOOLS wrapper that returns a `ZDOTC` result
through an output argument. `ZZDOTC` itself is defined in
`external/scalapack-2.2.3/TOOLS/zzdotc.f`, which `recipes/ptzblas.yaml`
does **not** list under `extra_symbol_dirs`, so:
- The migrator's PTZBLAS scan never sees `zzdotc` as a symbol, so the
  `CALL ZZDOTC(...)` line in the kind16 migration of `zvvdotc.f` isn't
  renamed to `CALL XDOTC(...)` the way it should be.
- `zzdotc.f` itself isn't staged or compiled into any migrated
  library, so the dangling call has nowhere to resolve.

**Local workaround.** `target_kind10/ptzblas_stubs.f90` and
`target_kind16/ptzblas_stubs.f90` provide thin `zzdotc` / `zzdotu`
stubs that forward to the target's canonical `ydotc`/`ydotu` and
`xdotc`/`xdotu` respectively. (For multifloats, the Z→ZZ double
substitution produces `zzzzdotc` naturally, so no stub is needed —
or at least none has been needed so far.)

**Fix upstream.** The cleanest fix is to add
`external/scalapack-2.2.3/TOOLS` to `recipes/ptzblas.yaml`'s
`extra_symbol_dirs`, and copy-and-compile `zzdotc.f` / `zzdotu.f`
into the staged `ptzblas/src/` tree so the migrated wrappers exist
alongside the renamed kernels. A narrower alternative is to patch
the PTZBLAS sources to call the direct BLAS function (`ZDOTC` →
`XDOTC`/`YDOTC`/`ZZDOTC` by renaming) rather than going through the
ScaLAPACK-TOOLS wrapper at all.

### 5. `pqsyrk` / `pxherk` produce numerically wrong results

**Symptom.** `pblas_test_pdsyrk` and `pblas_test_pzherk` fail with
max_rel_err around 1–5 % — i.e. only ~1–2 decimal digits of agreement
between the migrated PBLAS routine and serial quad-precision BLAS.
Every other Level-1/2/3 PBLAS routine in the suite agrees bit-exactly
at kind16.

**Reproduction.** `ctest -R '^pblas_test_p(d|z)(syrk|herk)'`. The
failure persists at both NP=4 (2×2 grid) and NP=1 (1×1 grid), so it's
not a distribution-layer issue — the migrated `pqsyrk_` / `pxherk_`
C entry points themselves are computing the wrong answer. Sample
output from the JSON reports:
```
pdsyrk kind16  n=32 k=24   max_rel_err=4.7e-02   digits=1.3 ✗
pdsyrk kind16  n=80 k=48   max_rel_err=2.4e-02   digits=1.6 ✗
pdsyrk kind16  n=128 k=100 max_rel_err=1.2e-02   digits=1.9 ✗
pzherk kind16  n=32 k=20   max_rel_err=3.5e-02   digits=1.5 ✗
pzherk kind16  n=64 k=40   max_rel_err=2.2e-02   digits=1.7 ✗
pzherk kind16  n=96 k=72   max_rel_err=1.2e-02   digits=1.9 ✗
```

**Status.** Not worked around — the test suite deliberately leaves
these as hard failures because that's what differential precision
tests are *for*. Likely culprits to investigate upstream: the
precision-family rename in `PB_Cpsyrk*` or the C↔Fortran scalar
passing for `ALPHA`/`BETA` in the herk/syrk code paths, or an
un-renamed call into a single-precision kernel inside the SUMMA
implementation. The error shrinks (4.7 % → 1.2 %) as `n` grows, which
hints at a fixed per-block contribution being wrong rather than a
uniform rounding issue.
