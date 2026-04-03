# BLAS Type Migration Recipe

This document describes how to migrate the Reference BLAS
(`external/lapack-3.12.1/BLAS/SRC/`) from standard precision
(`REAL`/`DOUBLE PRECISION`, `COMPLEX`/`COMPLEX*16`) to extended precision
(`KIND=10` or `KIND=16`).

## Source Overview

- **Location:** `external/lapack-3.12.1/BLAS/SRC/`
- **File count:** 155 Fortran source files (147 fixed-form `.f`, 8 free-form `.f90`)
- **Preprocessor:** None — pure Fortran throughout.

## File Classification

### Standard Precision-Prefixed Routines (Migratable)

These follow the S/D/C/Z prefix convention and are migrated by changing the
prefix and converting all precision-dependent constructs.

#### Single Real (S → E/Q): 33 files

```
sasum.f saxpy.f scopy.f sdot.f sgbmv.f sgemm.f sgemmtr.f sgemv.f sger.f
ssbmv.f sscal.f sspmv.f sspr.f sspr2.f sswap.f ssymm.f ssymv.f ssyr.f
ssyr2.f ssyr2k.f ssyrk.f stbmv.f stbsv.f stpmv.f stpsv.f strmm.f strmv.f
strsm.f strsv.f srot.f srotm.f srotmg.f snrm2.f90
```

#### Double Real (D → E/Q): 33 files

```
dasum.f daxpy.f dcopy.f ddot.f dgbmv.f dgemm.f dgemmtr.f dgemv.f dger.f
dsbmv.f dscal.f dspmv.f dspr.f dspr2.f dswap.f dsymm.f dsymv.f dsyr.f
dsyr2.f dsyr2k.f dsyrk.f dtbmv.f dtbsv.f dtpmv.f dtpsv.f dtrmm.f dtrmv.f
dtrsm.f dtrsv.f drot.f drotm.f drotmg.f dnrm2.f90
```

#### Single Complex (C → Y/X): 30 files

```
caxpy.f ccopy.f cdotc.f cdotu.f cgbmv.f cgemm.f cgemmtr.f cgemv.f cgerc.f
cgeru.f chbmv.f chemm.f chemv.f cher.f cher2.f cher2k.f cherk.f chpmv.f
chpr.f chpr2.f cscal.f cswap.f csymm.f csyr2k.f csyrk.f ctbmv.f ctbsv.f
ctpmv.f ctpsv.f ctrmm.f ctrmv.f ctrsm.f ctrsv.f crotg.f90
```

Note: `csrot.f` and `csscal.f` are cross-type routines (see below).

#### Double Complex (Z → Y/X): 30 files

```
zaxpy.f zcopy.f zdotc.f zdotu.f zgbmv.f zgemm.f zgemmtr.f zgemv.f zgerc.f
zgeru.f zhbmv.f zhemm.f zhemv.f zher.f zher2.f zher2k.f zherk.f zhpmv.f
zhpr.f zhpr2.f zscal.f zswap.f zsymm.f zsyr2k.f zsyrk.f ztbmv.f ztbsv.f
ztpmv.f ztpsv.f ztrmm.f ztrmv.f ztrsm.f ztrsv.f zrotg.f90
```

Note: `zdrot.f` and `zdscal.f` are cross-type routines (see below).

### Cross-Type / Mixed-Precision Routines (Special Handling)

These routines mix two precision types. They need careful analysis.

| File       | Description                  | Migration Notes                                |
|------------|------------------------------|------------------------------------------------|
| `dsdot.f`  | D result from S·S            | Mixed: `DOUBLE PRECISION` result, `REAL` args. Both types become the target kind, so rename to target prefix and convert all types. |
| `sdsdot.f` | S result with D accumulator  | Mixed: `REAL` result, local `DOUBLE PRECISION` accumulator. Both become target kind; accumulator becomes same type as result. |
| `csrot.f`  | Apply real rotation to C vec | Mixed: `COMPLEX` arrays, `REAL` rotation params. Both become target kind. |
| `csscal.f` | Scale C vector by S scalar   | Mixed: `COMPLEX` array, `REAL` scalar. Both become target kind. |
| `zdrot.f`  | Apply real rotation to Z vec | Mixed: `COMPLEX*16` arrays, `DOUBLE PRECISION` rotation params. Both become target kind. |
| `zdscal.f` | Scale Z vector by D scalar   | Mixed: `COMPLEX*16` array, `DOUBLE PRECISION` scalar. Both become target kind. |
| `dcabs1.f` | |Re(z)|+|Im(z)| for Z        | `DOUBLE PRECISION FUNCTION` with `COMPLEX*16` arg. Both become target kind. Uses `DBLE` and `DIMAG` intrinsics. |
| `scabs1.f` | |Re(z)|+|Im(z)| for C        | `REAL FUNCTION` with `COMPLEX` arg. Both become target kind. |

### Cross-Type Helper Routines (Return INTEGER)

These return INTEGER but take precision-specific array arguments.

| File        | Description                      | Migration Notes                        |
|-------------|----------------------------------|----------------------------------------|
| `icamax.f`  | Index of max |Re|+|Im| in C vec  | `COMPLEX` array arg → target kind. Calls `SCABS1`. |
| `izamax.f`  | Index of max |Re|+|Im| in Z vec  | `COMPLEX*16` array arg → target kind. Calls `DCABS1`. |
| `isamax.f`  | Index of max |val| in S vec      | `REAL` array arg → target kind. |
| `idamax.f`  | Index of max |val| in D vec      | `DOUBLE PRECISION` array arg → target kind. |

### Cross-Type Norm Routines (Return Real from Complex Input)

| File         | Description              | Migration Notes                        |
|--------------|--------------------------|----------------------------------------|
| `scnrm2.f90` | S norm of C vector      | `REAL` result, `COMPLEX` input → both become target kind. |
| `dznrm2.f90` | D norm of Z vector      | `DOUBLE PRECISION` result, `COMPLEX*16` input → both become target kind. Uses `real(wp)` and `complex(wp)` parameterized types. |
| `scasum.f`   | S sum of |Re|+|Im| of C | `REAL` result, `COMPLEX` input → both become target kind. |
| `dzasum.f`   | D sum of |Re|+|Im| of Z | `DOUBLE PRECISION` result, `COMPLEX*16` input → both become target kind. Calls `DCABS1`. |

### Precision-Independent Routines (No Migration Needed)

| File             | Description                  | Notes                                  |
|------------------|------------------------------|----------------------------------------|
| `lsame.f`        | Case-insensitive char compare| No floating-point types. Copy as-is.   |
| `xerbla.f`       | Error handler                | No floating-point types. Copy as-is.   |
| `xerbla_array.f` | Array-form error handler     | No floating-point types. Copy as-is.   |

### Fortran 90 Routines (Free-Form)

These 8 files use `.f90` free-form with parameterized types via
`integer, parameter :: wp = kind(1.d0)` and `real(wp)` / `complex(wp)`:

```
srotg.f90  drotg.f90  crotg.f90  zrotg.f90
snrm2.f90  dnrm2.f90  scnrm2.f90 dznrm2.f90
```

**Special migration considerations for `.f90` files:**
- The `wp` parameter definition (`kind(1.d0)` or `kind(1.e0)`) must be changed
  to the target kind value (e.g., `integer, parameter :: wp = 16`).
- Constants like `0.0_wp`, `1.0_wp` automatically follow `wp` — no literal
  conversion needed once `wp` is updated.
- Generic intrinsics (`abs`, `sqrt`, `sign`, `real`, `aimag`, `conjg`) are
  already used — no intrinsic conversion needed.
- Function/subroutine names still need prefix renaming.

---

## Transformation Details

### 1. Type Declarations

**Fixed-form files (`.f`):**

| Source Pattern         | → KIND=10          | → KIND=16          |
|------------------------|--------------------|--------------------|
| `REAL`                 | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `DOUBLE PRECISION`     | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `COMPLEX`              | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `COMPLEX*16`           | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |

Notes:
- `REAL` (6 → 13 chars for KIND=10, or 13 for KIND=16) can overflow column 72.
  Reformat into continuation lines if needed.
- `DOUBLE PRECISION` (16 chars) → `REAL(KIND=16)` (12 chars) is shorter —
  adjust trailing whitespace to preserve column alignment.
- `COMPLEX*16` (10 chars) → `COMPLEX(KIND=16)` (15 chars) may overflow.
- `COMPLEX` (7 chars) → `COMPLEX(KIND=16)` (15 chars) may overflow.
- Type keywords appear in declarations (column 7+), FUNCTION return types,
  PARAMETER statements, and EXTERNAL declarations with explicit types.
- Do NOT convert `INTEGER`, `LOGICAL`, `CHARACTER` — they are not
  precision-dependent.

**Free-form files (`.f90`):**

| Source Pattern              | → KIND=10                        | → KIND=16                        |
|-----------------------------|----------------------------------|----------------------------------|
| `integer, parameter :: wp = kind(1.d0)` | `integer, parameter :: wp = 10` | `integer, parameter :: wp = 16` |
| `integer, parameter :: wp = kind(1.e0)` | `integer, parameter :: wp = 10` | `integer, parameter :: wp = 16` |
| `DOUBLE PRECISION`          | `REAL(KIND=10)`                  | `REAL(KIND=16)`                  |
| `DOUBLE COMPLEX`            | `COMPLEX(KIND=10)`               | `COMPLEX(KIND=16)`               |

`real(wp)` and `complex(wp)` declarations do not need changing — they
automatically follow the `wp` parameter.

### 2. Subroutine/Function/Entry Names

Rename definitions, CALL sites, EXTERNAL declarations, and function references
in expressions.

**Prefix mapping (KIND=16):**

| Source | Target | Example          |
|--------|--------|------------------|
| S      | Q      | `SGEMM` → `QGEMM` |
| D      | Q      | `DGEMM` → `QGEMM` |
| C      | X      | `CGEMM` → `XGEMM` |
| Z      | X      | `ZGEMM` → `XGEMM` |

**Prefix mapping (KIND=10):**

| Source | Target | Example          |
|--------|--------|------------------|
| S      | E      | `SGEMM` → `EGEMM` |
| D      | E      | `DGEMM` → `EGEMM` |
| C      | Y      | `CGEMM` → `YGEMM` |
| Z      | Y      | `ZGEMM` → `YGEMM` |

**Where names appear:**
- `SUBROUTINE DGEMM(...)` → `SUBROUTINE QGEMM(...)`
- `COMPLEX*16 FUNCTION ZDOTC(...)` → `COMPLEX(KIND=16) FUNCTION XDOTC(...)`
- `DOUBLE PRECISION FUNCTION DSDOT(...)` → `REAL(KIND=16) FUNCTION QSDOT(...)`
  (Note: DSDOT is mixed-precision — see special handling)
- `CALL XERBLA('DGEMM ',INFO)` → `CALL XERBLA('QGEMM ',INFO)` (string literal)
- `EXTERNAL LSAME` — LSAME is type-independent, do NOT rename
- `EXTERNAL DCABS1` → `EXTERNAL QCABS1` (precision-specific external)
- `STEMP = STEMP + DCABS1(ZX(I))` → `STEMP = STEMP + QCABS1(ZX(I))`
- Comments: `*     End of DGEMM` → `*     End of QGEMM`
- Doc header: `*> \brief \b DGEMM` → `*> \brief \b QGEMM`

**Cross-type routine renaming:**

| Source     | KIND=16 | KIND=10 |
|------------|---------|---------|
| `DSDOT`    | `QSDOT`*| `ESDOT`*|
| `SDSDOT`   | `QSDOT`*| `ESDOT`*|
| `CSROT`    | `XSROT`*| `YSROT`*|
| `CSSCAL`   | `XSSCAL`*| `YSSCAL`*|
| `ZDROT`    | `XDROT`*| `YDROT`*|
| `ZDSCAL`   | `XDSCAL`*| `YDSCAL`*|
| `DCABS1`   | `QCABS1`*| `ECABS1`*|
| `SCABS1`   | `QCABS1`*| `ECABS1`*|
| `ICAMAX`   | `IXAMAX`*| `IYAMAX`*|
| `IZAMAX`   | `IXAMAX`*| `IYAMAX`*|
| `ISAMAX`   | `IQAMAX`*| `IEAMAX`*|
| `IDAMAX`   | `IQAMAX`*| `IEAMAX`*|
| `SCNRM2`   | `QXNRM2`*| `EYNRM2`*|
| `DZNRM2`   | `QXNRM2`*| `EYNRM2`*|
| `SCASUM`   | `QXASUM`*| `EYASUM`*|
| `DZASUM`   | `QXASUM`*| `EYASUM`*|

*The exact renaming of cross-type routines depends on the symbol database
classification. The prefix character that represents the primary precision
type is replaced. For routines with two precision characters (like `DZ`),
only the characters that are precision prefixes are changed. The symbol
database must be loaded to resolve these correctly.*

**Important: DO NOT rename these precision-independent routines:**
- `LSAME` — character comparison utility
- `XERBLA` — error reporting (the `X` is not a precision prefix)
- `XERBLA_ARRAY` — error reporting variant

### 3. Literal Constants

**Fixed-form (`.f`) patterns:**

| Source Literal        | → KIND=16         | → KIND=10         |
|-----------------------|-------------------|--------------------|
| `1.0D+0`              | `1.0E+0_16`       | `1.0E+0_10`        |
| `0.0D0`               | `0.0E0_16`        | `0.0E0_10`         |
| `1.0E+0`              | `1.0E+0_16`       | `1.0E+0_10`        |
| `0.0E+0`              | `0.0E+0_16`       | `0.0E+0_10`        |
| `(1.0D+0,0.0D+0)`     | `(1.0E+0_16,0.0E+0_16)` | `(1.0E+0_10,0.0E+0_10)` |
| `(0.0d0,0.0d0)`       | `(0.0E0_16,0.0E0_16)` | `(0.0E0_10,0.0E0_10)` |

Rules:
- Replace `D` exponent letter with `E` (D is DOUBLE PRECISION-specific).
- Append `_<kind>` suffix.
- Handle both uppercase and lowercase exponent letters (`D`, `d`, `E`, `e`).
- Complex literal pairs `(real, imag)` — each component is converted independently.
- `1.E0` (without decimal digit after `.`) is a valid Fortran literal — handle it.
- Literals appearing in comments are converted for documentation consistency.

**Free-form (`.f90`) patterns:**

Literals use `_wp` suffix (e.g., `0.0_wp`, `1.0_wp`) which automatically
follows the `wp` parameter. No literal conversion needed for these files
once `wp` is updated — only the `wp` parameter definition changes.

### 4. Intrinsic Functions

Type-specific intrinsics used in BLAS `.f` files (found via grep):

| Source Intrinsic | Replacement           | Files Using It          |
|------------------|-----------------------|-------------------------|
| `DCONJG(z)`      | `CONJG(z)`           | zgemm, zgemv, zdotc, zhemm, zher, zher2, zher2k, zherk, zhpmv, zhpr, zhpr2, ztbmv, ztbsv, ztpmv, ztpsv, ztrmm, ztrmv, ztrsm, ztrsv |
| `DBLE(x)`        | `REAL(x,KIND=k)`     | dcabs1, dsdot, sdsdot, zdscal, zgbmv |
| `DIMAG(z)`       | `AIMAG(z)`           | dcabs1, zgbmv, zgemv |
| `DCMPLX(x,y)`    | `CMPLX(x,y,KIND=k)` | zgerc, zhbmv, zgbmv, zgemv |
| `DABS(x)`        | `ABS(x)`             | drotmg (via INTRINSIC DABS) |

Also update the corresponding `INTRINSIC` declarations:
- `INTRINSIC DCONJG` → `INTRINSIC CONJG`
- `INTRINSIC DBLE` → remove (REAL is generic, no declaration needed; or `INTRINSIC REAL`)
- `INTRINSIC DIMAG` → `INTRINSIC AIMAG`
- `INTRINSIC DCMPLX` → remove (`CMPLX` is generic)
- `INTRINSIC DABS` → `INTRINSIC ABS`

Generic intrinsics already used in single/complex routines (`ABS`, `AIMAG`,
`CONJG`, `CMPLX`, `REAL`, `MAX`, `MIN`, `MOD`, `SQRT`, `SIGN`) need NO
conversion — they work for all KINDs.

### 5. XERBLA String Arguments

Every BLAS routine that validates its arguments calls:
```fortran
      CALL XERBLA('DGEMM ',INFO)
```
The string literal contains the routine name (padded to 6 or 7 chars with
trailing spaces). This must be updated to the new name:
```fortran
      CALL XERBLA('QGEMM ',INFO)
```

### 6. Comment Updates

BLAS files have extensive documentation headers. Routine names appear in:
- `*> \brief \b DGEMM` → `*> \brief \b QGEMM`
- `*> DGEMM performs one of the matrix-matrix operations`
- `*     End of DGEMM`
- `*> ALPHA is DOUBLE PRECISION` → `*> ALPHA is REAL(KIND=16)`
- Prototype blocks in the header documentation

While comments are not semantically critical, they serve as API documentation
and should be updated for consistency.

### 7. File Renaming

Output filenames follow the same prefix mapping:
- `dgemm.f` → `qgemm.f` (KIND=16) or `egemm.f` (KIND=10)
- `zgemm.f` → `xgemm.f` (KIND=16) or `ygemm.f` (KIND=10)
- `lsame.f` → `lsame.f` (no change — precision-independent)
- `dznrm2.f90` → `qxnrm2.f90` (KIND=16)

---

## Convergence Test Pairs

After migration, the following pairs should produce identical output
(modulo comment differences where the original S/D or C/Z text differs):

### Real Pairs (33 pairs)

```
sasum.f  ↔ dasum.f   → qasum.f
saxpy.f  ↔ daxpy.f   → qaxpy.f
scopy.f  ↔ dcopy.f   → qcopy.f
sdot.f   ↔ ddot.f    → qdot.f
sgbmv.f  ↔ dgbmv.f   → qgbmv.f
sgemm.f  ↔ dgemm.f   → qgemm.f
sgemmtr.f↔ dgemmtr.f → qgemmtr.f
sgemv.f  ↔ dgemv.f   → qgemv.f
sger.f   ↔ dger.f    → qger.f
ssbmv.f  ↔ dsbmv.f   → qsbmv.f
sscal.f  ↔ dscal.f   → qscal.f
sspmv.f  ↔ dspmv.f   → qspmv.f
sspr.f   ↔ dspr.f    → qspr.f
sspr2.f  ↔ dspr2.f   → qspr2.f
sswap.f  ↔ dswap.f   → qswap.f
ssymm.f  ↔ dsymm.f   → qsymm.f
ssymv.f  ↔ dsymv.f   → qsymv.f
ssyr.f   ↔ dsyr.f    → qsyr.f
ssyr2.f  ↔ dsyr2.f   → qsyr2.f
ssyr2k.f ↔ dsyr2k.f  → qsyr2k.f
ssyrk.f  ↔ dsyrk.f   → qsyrk.f
stbmv.f  ↔ dtbmv.f   → qtbmv.f
stbsv.f  ↔ dtbsv.f   → qtbsv.f
stpmv.f  ↔ dtpmv.f   → qtpmv.f
stpsv.f  ↔ dtpsv.f   → qtpsv.f
strmm.f  ↔ dtrmm.f   → qtrmm.f
strmv.f  ↔ dtrmv.f   → qtrmv.f
strsm.f  ↔ dtrsm.f   → qtrsm.f
strsv.f  ↔ dtrsv.f   → qtrsv.f
srot.f   ↔ drot.f    → qrot.f
srotm.f  ↔ drotm.f   → qrotm.f
srotmg.f ↔ drotmg.f  → qrotmg.f
snrm2.f90↔ dnrm2.f90 → qnrm2.f90
```

### Complex Pairs (30 pairs)

```
caxpy.f  ↔ zaxpy.f   → xaxpy.f
ccopy.f  ↔ zcopy.f   → xcopy.f
cdotc.f  ↔ zdotc.f   → xdotc.f
cdotu.f  ↔ zdotu.f   → xdotu.f
cgbmv.f  ↔ zgbmv.f   → xgbmv.f
cgemm.f  ↔ zgemm.f   → xgemm.f
...etc (all c* ↔ z* pairs)
crotg.f90↔ zrotg.f90 → xrotg.f90
```

### Cross-Type Convergence Pairs

```
scabs1.f ↔ dcabs1.f  → qcabs1.f
csrot.f  ↔ zdrot.f   → xdrot.f (or xsrot.f — verify via symbol DB)
csscal.f ↔ zdscal.f  → xdscal.f (or xsscal.f — verify via symbol DB)
scasum.f ↔ dzasum.f  → qxasum.f
scnrm2.f90↔dznrm2.f90→ qxnrm2.f90
icamax.f ↔ izamax.f  → ixamax.f
isamax.f ↔ idamax.f  → iqamax.f
```

### Routines Without Convergence Pairs

| Routine  | Reason                                    |
|----------|-------------------------------------------|
| `dsdot.f`  | Mixed S input / D output — no S equivalent |
| `sdsdot.f` | S result with D accumulator — unique       |
| `srotg.f90`| Has pair `drotg.f90` — listed above        |

---

## Execution Steps

### Step 1: Build Symbol Database

```bash
# From source (no pre-built library needed for BLAS):
fortran-migrator --kind 16 \
    external/lapack-3.12.1/BLAS/SRC \
    output/blas-q

# OR from compiled library:
fortran-migrator --kind 16 \
    --library build-lapack/lib/libblas.a \
    external/lapack-3.12.1/BLAS/SRC \
    output/blas-q
```

### Step 2: Migrate

The tool processes each file:
1. Load symbol database (from source scan or library)
2. For each `.f` / `.f90` file in source dir:
   a. Scan for type declarations → generate transformations
   b. Scan for routine names (definition, call, external) → transformations
   c. Scan for literal constants → transformations
   d. Scan for type-specific intrinsics → transformations
   e. Scan for XERBLA string arguments → transformations
   f. Validate no overlapping transformations
   g. Apply all transformations (back-to-front by offset)
   h. Write output with renamed filename

### Step 3: Convergence Check

```bash
fortran-migrator --kind 16 --convergence-check \
    external/lapack-3.12.1/BLAS/SRC \
    output/blas-q-convergence
```

This migrates from both S→Q and D→Q (and C→X and Z→X), then diffs each
convergence pair. Differences fall into categories:
- **Comment-only diffs:** harmless (e.g., "single precision" vs "double precision" in comments)
- **Algorithm diffs:** intentional (e.g., different blocking, thresholds)
- **Missed conversions:** bugs in the migrator

### Step 4: Compile Test

```bash
cd output/blas-q
gfortran -c -fdefault-real-16 *.f *.f90
# OR
gfortran -c *.f *.f90  # Should work if all types use explicit KIND=16
```

### Step 5: Build Libraries

The migrated output is split into two libraries:
- `libblas_common.a` — type-independent routines (LSAME, XERBLA, XERBLA_ARRAY)
- `libqblas.a` — precision-specific routines (Q/X/I-prefix)

This separation applies to all libraries (LAPACK, ScaLAPACK, etc.):
the common library is shared across all precision variants and linked once.

```bash
# Archive common (type-independent) objects
ar rcs libblas_common.a lsame.o xerbla.o xerbla_array.o

# Archive precision-specific objects (everything else)
ar rcs libqblas.a q*.o x*.o i*.o

# Link order: -lqblas -lblas_common
```

---

## Fixed-Form Column Considerations

BLAS `.f` files use fixed-form with strict column layout:
- Columns 1-5: labels
- Column 6: continuation (`+` is the convention in BLAS)
- Columns 7-72: statement
- Column 73+: ignored

Type replacements that change string length:
- `REAL` (4 chars) → `REAL(KIND=16)` (12 chars): **+8 chars** — may overflow
- `COMPLEX` (7 chars) → `COMPLEX(KIND=16)` (15 chars): **+8 chars** — may overflow
- `DOUBLE PRECISION` (16 chars) → `REAL(KIND=16)` (12 chars): **-4 chars** — shrinks
- `COMPLEX*16` (10 chars) → `COMPLEX(KIND=16)` (15 chars): **+5 chars** — may overflow

When a line exceeds column 72 after replacement:
1. First try absorbing extra whitespace in the line
2. If still too long, split into continuation lines using `+` in column 6

The continuation character used in BLAS is consistently `+`.

---

## Edge Cases Specific to BLAS

1. **`COMPLEX*16` vs `DOUBLE COMPLEX`:** BLAS uses `COMPLEX*16` (non-standard
   but common) rather than `DOUBLE COMPLEX`. Both map to `COMPLEX(KIND=16)`.

2. **Complex literal pairs:** `(1.0D+0,0.0D+0)` — each real/imaginary component
   must be converted independently. The parentheses and comma are preserved.

3. **`DCONJG` in `INTRINSIC` list:** `INTRINSIC DCONJG,MAX` — the `DCONJG`
   must be replaced with `CONJG` within the comma-separated list.

4. **Function return types:** `COMPLEX*16 FUNCTION ZDOTC(...)` — both the
   return type (`COMPLEX*16`) and function name (`ZDOTC`) need conversion.

5. **`DROTMG`/`SROTMG`:** Contains extensive literal constants in comments
   (e.g., `1.E0`, `-1.E0`) and in code. All must be converted.

6. **Free-form `wp` parameter:** The `.f90` files define
   `integer, parameter :: wp = kind(1.d0)`. This single line controls all
   precision — changing it to `wp = 16` is sufficient for types and literals.
   Only routine names still need explicit renaming.

7. **`srotg.f90`:** Uses `kind(1.e0)` for single precision.
   `drotg.f90` uses `kind(1.d0)` for double precision. Both → `16` or `10`.

8. **LSAME, XERBLA:** These appear in `EXTERNAL` declarations in many files.
   They must NOT be renamed. The symbol database handles this by classifying
   them as `PrecisionKind::None`.
