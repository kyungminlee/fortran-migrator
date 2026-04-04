# fortran-migrator — Progress & Notes

## Branch: `claude/test-blacs-migration-i8XlV`

### What This Branch Does

Adds support for migrating BLACS (C library) and ScaLAPACK (Fortran) alongside
the existing BLAS/LAPACK migration, including proper MPI datatype handling,
Bdef.h patching, LAPACK constant module support, and intrinsic function fixes
for extended/quad precision.

---

## Pipeline Status (all tested on this branch)

| Library | KIND=16 | KIND=10 | Notes |
|---------|---------|---------|-------|
| BLAS    | migrate ✓ verify ✓ build ✓ | migrate ✓ verify ✓ build ✓ | Clean |
| LAPACK  | migrate ✓ verify ⚠ build ✓ | migrate ✓ verify ⚠ build ✓ | 6 benign verify warnings (see below) |
| BLACS   | migrate ✓ build ✓ | migrate ✓ build ✓ | No verify step for C yet |
| ScaLAPACK | migrate ✓ verify ⚠ build ✓ | migrate ✓ verify ⚠ build ✓ | 1 benign verify warning (see below) |

### LAPACK Verify Warnings (benign — do not need fixing)

1. **`la_constants.f90` residuals (2 warnings)**: Lines defining `sp = kind(1.e0)`
   and `dp = kind(1.d0)` — this is the *original* module that is copied unchanged
   (it's in `copy_files`). The migrated code uses `LA_CONSTANTS_EP` instead.

2. **`xlarf1f.f` column overflow (4 warnings)**: Lines 234, 275, 291, 294 exceed
   72 chars. These are on *comment lines* where the original `dlarf1f.f` was
   already close to the limit and the `d→q` prefix change pushed them over.
   The Fortran compiler ignores comment overflow.

### ScaLAPACK Verify Warning (benign — do not need fixing)

1. **`piparmq.f` D-exponent residual (1 warning)**: Line 263 has `335.0D+0`
   and `-0.44D+0` — this is a precision-independent file copied unchanged
   from the original source. The D-exponents are in the original code.

---

## Key Architecture Decisions

### LA_CONSTANTS_EP Module (`external/lapack-3.12.1/SRC/la_constants_ep.F90`)

LAPACK's `la_constants.f90` defines single-precision (S-prefix) and
double-precision (D-prefix) scaling constants. Rather than modifying it, we
created a **companion module** `LA_CONSTANTS_EP` with:
- Quad-precision (Q-prefix) constants — always available
- Extended-precision (E-prefix) constants — guarded by `#ifdef HAVE_REAL10`

The migrator rewrites `USE LA_CONSTANTS` → `USE LA_CONSTANTS_EP` and renames
all imported symbols (e.g., `wp=>dp` → `wp=>qp`, `DZERO` → `QZERO`).

### LA_XISNAN_EP Module (`external/lapack-3.12.1/SRC/la_xisnan_ep.F90`)

Same pattern: companion to `LA_XISNAN` adding `QISNAN` and `EISNAN` overloads
for the generic `LA_ISNAN` interface. Required because some LAPACK routines
(e.g., `qlassq.f90`) call `LA_ISNAN` with quad-precision arguments.

### USE Statement Rewriting (`rewrite_la_constants_use()` in `fortran_migrator.py`)

The symbol scanner only finds routine/function names, NOT module parameters.
So constants like `SP`, `DP`, `SZERO`, `DZERO` were never in the rename map.
We added a **hardcoded constant list** (`_LA_CONST_SP`, `_LA_CONST_DP`, etc.)
and a dedicated rewriting pass that:
1. Changes the module name (`LA_CONSTANTS` → `LA_CONSTANTS_EP`)
2. Renames all imported symbols using the prefix map
3. Handles multi-line USE statements (continuation with `&`)
4. Also rewrites `USE LA_XISNAN` → `USE LA_XISNAN_EP`

This runs as a **pre-processing step** before the line-by-line migration in
both `migrate_free_form()` and `_migrate_free_form_flang()`.

### ScaLAPACK Intrinsic Fixes (in `intrinsics.py` and `fortran_migrator.py`)

Three issues surfaced during ScaLAPACK testing that required migrator fixes:

1. **IFIX/IDINT intrinsics**: `IFIX` (REAL(4)→INTEGER) and `IDINT`
   (REAL(8)→INTEGER) are precision-specific integer conversion intrinsics
   used in `psgeqpf.f`, `pcgeqpf.f`, `pdgeqpf.f`, `pzgeqpf.f`. Both now
   map to the generic `INT`.

2. **Duplicate INTRINSIC declarations**: When `DBLE` → `REAL` replacement
   creates a duplicate in an INTRINSIC declaration (e.g., the original has
   both `DBLE` and `REAL`), gfortran rejects it. The `replace_intrinsic_decls()`
   function now deduplicates names while preserving trailing commas needed
   for fixed-form continuation lines.

3. **DBLE(variable) wrongly skipped**: The `replace_intrinsic_calls()` function
   had a heuristic to skip `REAL(wp)` type declarations (single identifier
   argument). This incorrectly skipped `DBLE(ALPHA)` function calls. Fixed by
   only applying the single-identifier skip for names that can be Fortran type
   specifiers (`REAL`, `CMPLX`, `COMPLEX`), not for `DBLE`, `SNGL`, etc.

### HAVE_REAL10 Macro

KIND=10 (80-bit extended) is not universally supported:
- **gfortran on x86**: supports it (`selected_real_kind(18, 4931)` returns 10)
- **Intel Fortran (ifx/ifort)**: does NOT support it
- **flang**: does NOT support it

The generated CMakeLists.txt uses `CheckFortranSourceCompiles` to auto-detect
support and defines `HAVE_REAL10` if available. The `.F90` files use
`#ifdef HAVE_REAL10` guards (note: uppercase `.F90` enables preprocessing).

### BLACS MPI Type Mapping

| KIND | Real MPI Type | Complex MPI Type | C Real Type | C Complex Type |
|------|--------------|-----------------|-------------|---------------|
| 16 | `MPI_REAL16` | `MPI_COMPLEX32` | `__float128` | `struct {__float128 r, i;}` |
| 10 | `MPI_LONG_DOUBLE` | `MPI_C_LONG_DOUBLE_COMPLEX` | `long double` | `struct {long double r, i;}` |

- `MPI_REAL16` / `MPI_COMPLEX32` are **not in the MPI standard** but are
  supported by OpenMPI and MPICH when built with `__float128` support.
- A `CheckMpiReal16.cmake` module is generated for KIND=16 builds that
  `FATAL_ERROR`s at configure time if the MPI implementation doesn't provide them.
- For KIND=10, `MPI_LONG_DOUBLE` is standard MPI and maps to `long double`
  (80-bit extended on x86).

### BLACS Bdef.h Patching

The C migrator post-processes `Bdef.h` in the output directory to add:
- `typedef __float128 QREAL` / `typedef long double EREAL`
- `typedef struct {QREAL r, i;} XCOMPLEX` / `YCOMPLEX`
- Function prototypes for `BI_{rp}mvcopy` / `BI_{rp}vmcopy`
- Complex copy macros (`BI_xmvcopy` delegates to `BI_qmvcopy`)
- Fortran name mangling `#define`s in NOCHANGE and UPCASE sections

The FCISF2C (f2c) section is NOT patched — it's rarely used and would need
double-underscore variants.

---

## Caveats & Known Limitations

### 1. BLACS has no `cmd_verify` support
The verify step (`pyengine verify`) only handles Fortran files. There is no
equivalent check for C files. Manual inspection is needed for BLACS.

### 2. BLACS has no `cmd_build` support
The `pyengine build` command generates a Fortran-only CMakeLists.txt. For
BLACS, you must provide your own CMakeLists.txt. The test builds used:
```cmake
find_package(MPI REQUIRED COMPONENTS C)
add_compile_definitions(Add_)
file(GLOB BLACS_SOURCES src/*.c)
add_library(blacs STATIC ${BLACS_SOURCES})
target_link_libraries(blacs PUBLIC MPI::MPI_C)
```

### 3. LAPACK `copy_files` must include EP modules
The recipe's `copy_files` list must include both original and EP modules:
```yaml
copy_files:
  - LA_CONSTANTS
  - LA_CONSTANTS_EP
  - LA_XISNAN
  - LA_XISNAN_EP
```
If you forget the EP entries, the migrated code will import from a module
that doesn't exist in the output directory.

### 4. `__float128` compiler support
- GCC: Full support for `__float128` in C
- Clang: Partial — works on x86-64 Linux but not all platforms
- The BLACS C code does arithmetic on `__float128` values in reduction
  helpers (`BI_qvvamx.c`, `BI_qvvsum.c`, etc.), so it's not just opaque
  byte passing — the compiler must support `__float128` math.

### 5. ScaLAPACK `piparmq.f` D-exponent residual
The copied (precision-independent) file `piparmq.f` contains `335.0D+0`
and `-0.44D+0` literals. This is in the original source and triggers a
benign verify warning. The file is not migrated, just copied.

### 6. Fixed-form Fortran column limits
The migrator does not wrap or reformat lines that exceed 72 columns after
prefix substitution. This occasionally produces column overflow in comments
(harmless) but could in theory affect code lines too. The `xlarf1f.f`
warnings are the only known instances and are all on comment lines.

### 7. Fortran name mangling in Bdef.h
Only the `NOCHANGE` and `UPCASE` sections are patched with q/x (or e/y)
routine defines. The `ADD_` section needs no defines (underscore suffix is
the default). The `FCISF2C` section (double underscores for f2c) is NOT
patched — add it if you need f2c support.

---

## Environment Tested

- **OS**: Ubuntu 24.04 (x86-64)
- **gfortran**: 13.3.0
- **gcc**: 13.3.0
- **OpenMPI**: 4.1.6 (MPI 3.1) — supports MPI_REAL16
- **CMake**: 3.28.3
- **flang-new**: 18.1.3 (LLVM) — available but only used for BLAS/LAPACK,
  not BLACS (BLACS is pure C)

---

## Commit History (this branch)

```
765b65c Fix BLACS C migration: rename Fortran-callable functions and BI_*MPI_* helpers
12b95f0 Add la_constants_ep.F90 for extended/quad precision LAPACK constants
1c52897 Rewrite USE LA_CONSTANTS/LA_XISNAN to EP modules in migrated code
44c6c23 Use standard MPI types (MPI_REAL16/MPI_COMPLEX32) for quad-precision BLACS
358e41e Patch Bdef.h with extended-precision types, macros, and name mangling
e69dd01 Add .agent/PROGRESS.md with project status and continuation notes
3387f6d Fix ScaLAPACK migration: add IFIX/IDINT intrinsics, dedup INTRINSIC decls, fix DBLE(var) skip
```

---

## Likely Next Steps

1. **Integrated build**: Build all four libraries together (BLAS → LAPACK →
   BLACS → ScaLAPACK) and link them into a single test program.
2. **C migrator `cmd_build` support**: Generate a proper CMakeLists.txt for
   C libraries (with MPI, compiler flags, etc.) from `pyengine build`.
3. **C migrator `cmd_verify`**: Add verification for C files (check for
   residual `double`/`MPI_DOUBLE` in q/x files).
4. **Test with Intel Fortran / NVIDIA HPC SDK**: Verify KIND=16 works with
   non-GCC compilers (KIND=10 will be skipped via HAVE_REAL10 guard).
