# Developer Guide

## Additional Documentation

*   [Usage Guide](USAGE.md) — How to run the `pyengine` CLI.
*   [Migration Recipes](RECIPES.md) — How to write library recipes.
*   [Architecture](ARCHITECTURE.md) — Technical overview and component details.
*   [Convergence Testing](CONVERGENCE_TESTING.md) — The dual-origin verification strategy.

## Project Goal

**fortran-migrator** converts BLAS, LAPACK, BLACS, ScaLAPACK, and similar
numerical Fortran libraries from standard precision (`REAL`/`DOUBLE PRECISION`,
`COMPLEX`/`DOUBLE COMPLEX`) to extended precision (`KIND=10` for 80-bit,
`KIND=16` for 128-bit). This involves:

1. Rewriting type declarations
2. Renaming subroutine/function names (both definitions and call sites)
3. Updating literal constants
4. Converting type-specific intrinsic calls
5. Updating `EXTERNAL` declarations
6. Renaming output file names

The tool must preserve all preprocessor directives, comments, and formatting.

---

## Approach: Hybrid (Flang Parse Tree + Source-Level Rewriting)

### Why Not Pure String Manipulation?

Fixed-form Fortran has complex rules (columns 1-5 labels, column 6 continuation,
columns 7-72 statements). A type declaration can span multiple continuation
lines, appear inside preprocessor-guarded blocks, or use non-standard syntax
like `REAL*8`, `COMPLEX*16`, `DOUBLE PRECISION`. Simple regex is brittle and
error-prone — it risks false positives in comments, string literals, or
identifiers that happen to contain type keywords.

### Why Not Pure Flang Unparse?

Flang's `flang-new -fdebug-unparse` can round-trip source code, but it
normalizes formatting, rewrites continuation lines, and may alter whitespace
in ways that make the output differ from the original. For libraries like
LAPACK that have very precise fixed-form formatting and extensive comment
blocks, this is unacceptable. Preprocessing directives can also be lost or
rearranged.

### The Hybrid Strategy

1. **Parse with Flang's parser library** (`flang/lib/Parser`) to produce a
   parse tree with full source provenance (exact byte offsets in original source).
2. **Walk the parse tree** to identify all locations that need transformation:
   type declarations, function/subroutine names, `EXTERNAL` statements, literal
   constants, intrinsic function calls, `CALL` statements, and string literals
   containing routine names (e.g., `CALL XERBLA('DGEMM ', INFO)`).
3. **Apply transformations as source-level byte-range replacements**, preserving
   everything else (comments, preprocessor directives, whitespace, column
   alignment) exactly as-is.
4. **Handle column-width constraints** for fixed-form: when a replacement string
   is longer than the original (e.g., `DOUBLE PRECISION` → `REAL(KIND=16)` is
   shorter, but `REAL` → `REAL(KIND=16)` is longer), adjust spacing or split
   into continuation lines as needed while staying within columns 7-72.

This gives us the accuracy of an AST-based approach with the preservation
guarantees of text-level patching.

### Alternative: Lightweight Custom Scanner (Fallback)

If linking against LLVM/Flang proves too heavyweight for distribution, a
fallback approach is a purpose-built Fortran scanner in C++ that understands
just enough of Fortran's lexical structure to identify:
- Type declaration statements (at the start of a declaration section)
- `SUBROUTINE` / `FUNCTION` / `ENTRY` names
- `CALL` targets
- `EXTERNAL` listed names
- Literal constants (`1.0D+0`, `1.0E+0`)
- Intrinsic function names in `INTRINSIC` statements and call expressions

This scanner would not need to be a full parser — it only needs to distinguish
statements from comments and string literals, handle continuation lines, and
recognize the specific syntactic patterns we care about. It can be significantly
simpler than a full Fortran parser.

---

## What Gets Transformed

### 1. Type Declarations

Source precision types and their target mappings:

| Source Type            | KIND=10 Target     | KIND=16 Target     |
|------------------------|--------------------|--------------------|
| `REAL`                 | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `REAL*4`               | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `REAL*8`               | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `DOUBLE PRECISION`     | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `REAL(KIND=8)`         | `REAL(KIND=10)`    | `REAL(KIND=16)`    |
| `COMPLEX`              | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `COMPLEX*8`            | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `COMPLEX*16`           | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `DOUBLE COMPLEX`       | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |
| `COMPLEX(KIND=8)`      | `COMPLEX(KIND=10)` | `COMPLEX(KIND=16)` |

Note: the source files may use *either* single-precision (`REAL`, `COMPLEX`)
or double-precision (`DOUBLE PRECISION`, `DOUBLE COMPLEX`, `REAL*8`,
`COMPLEX*16`) depending on the routine variant (`s*`/`c*` vs `d*`/`z*`).
All floating-point types get mapped to the target kind.

### 2. Subroutine/Function Name Prefixes

Standard BLAS/LAPACK naming convention uses a single-character prefix:

| Source Prefix | Meaning           | KIND=10 Prefix | KIND=16 Prefix |
|---------------|-------------------|----------------|----------------|
| `S`           | Single real       | `E`            | `Q`            |
| `D`           | Double real       | `E`            | `Q`            |
| `C`           | Single complex    | `Y`            | `X`            |
| `Z`           | Double complex    | `Y`            | `X`            |

Names must be updated in:
- `SUBROUTINE` / `FUNCTION` / `ENTRY` definition lines
- `CALL` statements
- `EXTERNAL` declarations
- Function references in expressions (e.g., `TEMP = DDOT(N, X, 1, Y, 1)`)
- String literals passed to `XERBLA` (e.g., `'DGEMM '` → `'QGEMM '`)

For ScaLAPACK, the prefix is `P` + precision character (e.g., `PDGESV` →
`PQGESV`). The `P` prefix is preserved; only the precision character changes.

### 3. Literal Constants

| Source Literal  | KIND=10 Replacement    | KIND=16 Replacement    |
|-----------------|------------------------|------------------------|
| `1.0E+0`        | `1.0E+0_10`*           | `1.0E+0_16`*           |
| `1.0D+0`        | `1.0E+0_10`*           | `1.0E+0_16`*           |
| `0.0D0`         | `0.0E0_10`*            | `0.0E0_16`*            |

\* The exact syntax for extended-precision literals in Fortran is
`value_kind` (e.g., `1.0_16`). The `D` exponent letter is specific to
`DOUBLE PRECISION` and must be replaced with `E` when using explicit `KIND`.
Alternatively, for `KIND=16`, the `Q` exponent letter is a common extension
(`1.0Q+0`), but this is non-standard and compiler-dependent. Prefer the
standard `_kind` suffix form.

### 4. Intrinsic Functions

Type-specific intrinsics must be converted to their generic or
extended-precision equivalents:

| Source Intrinsic | Generic Equivalent | Notes |
|------------------|-------------------|-------|
| `DBLE(x)`        | `REAL(x, KIND=k)` | Type conversion |
| `REAL(x)`        | `REAL(x, KIND=k)` | Add KIND argument |
| `DCMPLX(x,y)`   | `CMPLX(x, y, KIND=k)` | |
| `DCONJG(z)`      | `CONJG(z)` | Generic works for all kinds |
| `DIMAG(z)`       | `AIMAG(z)` | Generic works for all kinds |
| `DABS(x)`        | `ABS(x)` | Generic |
| `DSQRT(x)`       | `SQRT(x)` | Generic |
| `DEXP(x)`        | `EXP(x)` | Generic |
| `DLOG(x)`        | `LOG(x)` | Generic |
| `DSIN(x)` / `DCOS(x)` | `SIN(x)` / `COS(x)` | Generic |
| `DSIGN(x,y)`     | `SIGN(x,y)` | Generic |
| `DMAX1(x,y)`     | `MAX(x,y)` | Generic |
| `DMIN1(x,y)`     | `MIN(x,y)` | Generic |
| `DNINT(x)`       | `ANINT(x)` | Generic |
| `IDNINT(x)`      | `NINT(x)` | Generic |

### 5. Machine-Parameter Routines

`DLAMCH` / `SLAMCH` provide machine epsilon, overflow/underflow thresholds, etc.
These routines must be regenerated for the target precision — they cannot simply
be renamed because the constants they return are precision-dependent. The
migrated version (`QLAMCH` / `ELAMCH` / `XLAMCH` / `YLAMCH`) needs to compute
and return the correct machine parameters for the target `KIND`.

Similarly, `DLAMC3` and other helper routines used by `DLAMCH` need attention.

### 6. File Names

Output files are renamed to reflect the new prefix:
- `dgemm.f` → `qgemm.f`
- `zgesv.f` → `xgesv.f`
- `pdgesv.f` → `pqgesv.f`

---

## Symbol Database

The tool uses a **symbol database** built from compiled libraries to drive
renaming decisions. This ensures only functions that actually exist in the
library are renamed, avoiding false positives on identifiers that happen to
start with `d`/`s`/`c`/`z`.

### Building the Symbol Database

1. **Extract symbols** from a built static library (`.a`) or shared library
   (`.so`/`.dylib`) using `nm`:
   ```
   nm --defined-only liblapack.a | grep ' [TtWw] ' | awk '{print $3}'
   ```
2. **Normalize** symbol names: strip leading underscore (macOS convention),
   strip trailing underscore (Fortran convention), convert to uppercase.
3. **Classify** each symbol by its precision prefix to build a mapping:
   `DGEMM` maps to `{base: "GEMM", prefix: "D", type: "real"}`.
4. **Generate rename map**: for each symbol, compute the target name using the
   prefix mapping table.

### Handling Libraries Without Pre-Built Binaries

For source-only distributions, the tool can scan the source directory for
`SUBROUTINE` and `FUNCTION` definitions to build the symbol list. This is
less reliable (conditional compilation may hide some definitions) but works
for the common case.

---

## Target Libraries

### BLAS (Basic Linear Algebra Subprograms)
- **Source format**: Fixed-form (`.f`)
- **Naming**: Single-character prefix (`DGEMM`, `SGEMV`, `ZHER2K`)
- **Location**: `external/lapack-3.12.1/BLAS/SRC/`
- **Notes**: Pure Fortran, no preprocessor usage in source files.
  Uses `COMPLEX*16` rather than `DOUBLE COMPLEX` in some routines.

### LAPACK (Linear Algebra PACKage)
- **Source format**: Fixed-form (`.f`)
- **Naming**: Single-character prefix (`DGESV`, `ZHEEV`, `SGEQRF`)
- **Location**: `external/lapack-3.12.1/SRC/`
- **Notes**: Large library (~2000 source files). Pure Fortran,
  no preprocessor. Heavy use of `EXTERNAL` declarations and cross-routine calls.

### BLACS (Basic Linear Algebra Communication Subprograms)
- **Source format**: C (`.c`) with headers
- **Location**: `external/scalapack-2.2.3/BLACS/SRC/`
- **Notes**: BLACS is implemented entirely in C. Type-specific routines use
  `d`/`s`/`c`/`z` in C function names. Migration is done via **templated
  cloning** (clone `d*` → `q*`, `z*` → `x*` with mechanical text
  substitution) rather than AST-based transformation. The C code is
  structurally simpler than Fortran — no column constraints, no multi-word
  type keywords, and type-specific files are near-identical clones of each
  other. The header `Bdef.h` must be extended to add `QCOMPLEX`/`QREAL`
  typedefs, copy macros, and Fortran naming `#define` entries. Custom MPI
  datatypes must be registered for the extended-precision types. A clang
  parser is not needed. See `recipes/blacs/README.md` for the full recipe.

### ScaLAPACK (Scalable LAPACK)
- **Source format**: Fixed-form Fortran (`.f`) and some C (`.c`)
- **Naming**: `P` + precision prefix (`PDGESV`, `PZHEEV`)
- **Location**: `external/scalapack-2.2.3/SRC/`
- **Notes**: Follows same patterns as LAPACK but with distributed-memory
  descriptors. Calls BLACS and PBLAS routines. Some C source files exist
  alongside Fortran.

### PBLAS (Parallel BLAS)
- **Source format**: Mixed C and Fortran
- **Location**: `external/scalapack-2.2.3/PBLAS/SRC/`
- **Notes**: C wrapper layer around BLAS with MPI. Uses preprocessor macros
  extensively.

### MUMPS (future)
- **Naming**: Prefix + `MUMPS` (`DMUMPS`, `ZMUMPS`)
- **Notes**: Same `S`/`D`/`C`/`Z` prefix convention. Free-form Fortran (`.f90`).
  Uses preprocessor directives extensively. Good test case for free-form support.

---

## Fortran Format Considerations

### Fixed-Form (`.f`, `.for`)
- **Columns 1-5**: Statement label (numeric)
- **Column 6**: Continuation character (any non-blank, non-zero character; typically `+`, `$`, `&`, or `*`)
- **Columns 7-72**: Statement body
- **Column 73+**: Ignored (historically sequence numbers)
- **Comments**: `C`, `c`, `*`, or `!` in column 1
- Continuation lines in BLAS/LAPACK typically use `+` or `$` in column 6

When a replacement makes a line exceed column 72, the tool must:
1. Try to reduce whitespace within the statement
2. If still too long, split into continuation lines
3. Preserve the continuation character style used in the file

### Free-Form (`.f90`, `.f95`, `.F90`)
- No column restrictions (but typically 132 characters max)
- `&` at end of line for continuation; `&` at start of next line optional
- `!` for inline comments
- Case insensitive (like all Fortran)

### Preprocessor Directives
- Lines starting with `#` (after optional whitespace in free-form)
- Must be preserved verbatim — never modified
- Code inside `#ifdef` blocks should still be transformed
- `#include` files may need separate processing

---

## Architecture

### Pipeline

```
Input Files ──→ [Symbol DB] ──→ [Scanner/Parser] ──→ [Transformation Plan] ──→ [Rewriter] ──→ Output Files
                    │                  │                       │
                 nm / source        Flang or             List of
                 scan              custom scanner      (offset, len, replacement)
```

### Components

1. **SymbolDatabase**
   - Input: library archive (`.a`) or source directory
   - Output: set of known routine names with prefix classification
   - Builds rename map from source prefix to target prefix

2. **SourceScanner**
   - Input: single Fortran source file
   - Output: list of `Transformation` records `{offset, length, replacement_text}`
   - Uses Flang parse tree (or custom scanner) to identify all sites
   - Must handle both fixed-form and free-form

3. **TransformationEngine**
   - Validates that transformations don't overlap
   - Sorts by offset (descending) for safe in-place replacement
   - Handles column-width adjustments for fixed-form

4. **FileRewriter**
   - Applies transformations to produce output text
   - Writes output file with renamed filename

5. **Driver / CLI**
   - Accepts: source directory, output directory, target KIND, library archive
   - Orchestrates the pipeline
   - Supports dry-run mode (show what would change)
   - Parallel processing of independent files

---

## Build Dependencies

| Dependency     | Purpose | Required? |
|----------------|---------|-----------|
| CMake ≥ 3.27   | Build system | Yes |
| C++20 compiler | Implementation language | Yes |
| LLVM/Flang     | Fortran parser library (`libflangParser`, `libflangCommon`, `libflangSupport`) | Yes (for Flang-based approach) |
| `nm` / `objdump` | Symbol extraction from compiled libraries | Optional (system tool) |

### Linking Against Flang

Flang's parser is part of the LLVM monorepo. The relevant libraries are:

- `flangParser` — tokenizer, prescan, parse tree
- `flangCommon` — Fortran kind parameters, character utilities
- `flangSupport` — source provenance, encoding
- `flangEvaluate` — expression semantics (may be useful for constant folding)

These can be found via CMake's `find_package(Flang)` or by pointing to an
LLVM build directory.

---

## Edge Cases and Pitfalls

1. **XERBLA string arguments**: Routine names appear as string literals in
   `CALL XERBLA('DGEMM ', INFO)`. These must be renamed too.

2. **Mixed-precision routines**: Some routines call across precisions (e.g.,
   `DSDOT` calls both single and double). These need careful analysis — the
   routine itself may not have a direct quad equivalent.

3. **Machine constants**: `DLAMCH` returns precision-specific constants.
   The migrated version must recompute these for the target precision.

4. **SAVE'd variables and DATA statements**: Type-specific constants in `DATA`
   statements must be updated.

5. **EQUIVALENCE statements**: These create type aliasing that may break if
   type sizes change. Flag these for manual review.

6. **Fortran `INCLUDE`**: Some files use `INCLUDE 'filename'` which is not a
   preprocessor directive but a Fortran statement. Included files may contain
   type declarations that need migration.

7. **Comment-only name references**: Names in comments (especially the header
   documentation blocks) should ideally be updated for consistency but are not
   semantically critical. The tool should update them since LAPACK uses these
   comments as API documentation.

8. **`PARAMETER` constants**: e.g., `PARAMETER (ONE=1.0D+0, ZERO=0.0D+0)` —
   both the literal format and the variable type must be consistent.

9. **C source files (BLACS, CBLAS, LAPACKE)**: BLACS is implemented
   entirely in C and uses a template-clone pattern — type-specific files
   (`dgebr2d_.c`, `sgebr2d_.c`, etc.) are near-identical and differ only
   in C type names (`double`→`float`) and MPI datatype constants. Migration
   is done by cloning `d*`→`q*` / `z*`→`x*` with mechanical text
   substitution on a fixed set of tokens (`double`→`QREAL`,
   `MPI_DOUBLE`→`MPI_QREAL`, etc.). No clang parser is needed because C
   types are unambiguous single tokens, there are no column constraints,
   and no floating-point literal format issues (C has no `D` exponent).
   New extended-precision type definitions (`QCOMPLEX`, `QREAL`) and
   custom MPI datatypes are added to the `Bdef.h` header.
   CBLAS and LAPACKE are C wrappers around Fortran routines and follow
   a similar clone pattern.

10. **Non-prefixed helper routines**: Some LAPACK internal routines don't
    follow the prefix convention (e.g., `LSAME`, `XERBLA`, `ILAENV`).
    These are type-independent and should *not* be renamed. The symbol database
    handles this by only renaming symbols that have confirmed precision variants.

---

## Testing Strategy

### 1. Round-Trip Test
Parse → identify → apply null transformation → verify output is byte-identical
to input.

### 2. Known-Good Conversions
Manually create expected output for a small set of files (`dgemm.f` → `qgemm.f`,
`zgesv.f` → `xgesv.f`) and verify the tool produces identical results.

### 3. Dual-Origin Convergence Test

Every BLAS/LAPACK routine exists in both single- and double-precision variants
(e.g., `sgemm.f` and `dgemm.f`). When both are migrated to the same target
precision, they should produce **identical** output (modulo the source origin
comment, if any). This is a powerful validation technique:

```
sgemm.f ──migrate(KIND=16)──→ qgemm_from_s.f ─┐
                                                ├── diff
dgemm.f ──migrate(KIND=16)──→ qgemm_from_d.f ─┘
```

Similarly for complex:
```
cgemm.f ──migrate(KIND=16)──→ xgemm_from_c.f ─┐
                                                ├── diff
zgemm.f ──migrate(KIND=16)──→ xgemm_from_d.f ─┘
```

**What agreement tells us**: the migrator correctly normalizes all
precision-dependent constructs — type declarations, literal constants,
intrinsic calls, and routine names — regardless of the starting precision.

**What differences reveal**:
- **Algorithm differences**: Some routines have genuinely different
  implementations between single and double precision (e.g., different
  blocking factors, different convergence thresholds, extra iteration
  steps). These are intentional and the diff highlights them for review.
- **Hardcoded constants**: A literal `1.0E-6` in a single-precision routine
  vs `1.0D-12` in the double-precision version indicates a precision-dependent
  tolerance that needs manual attention for the target precision.
- **Missing intrinsic conversions**: If `SNRM2` became `QNRM2` but `DNRM2`
  also became `QNRM2` via a different path, the result should match. A
  mismatch means one path missed a conversion.
- **Mixed-precision routines**: Routines like `DSDOT` (double-precision result
  from single-precision inputs) or `ZCGESV` (mixed single/double complex
  iterative refinement) have no single-precision counterpart. These will
  naturally be flagged as having no convergence pair.
- **Comment/documentation drift**: Different wording in comment blocks
  between `s*` and `d*` variants is harmless but shows up in the diff.
  These can be filtered out or accepted.

**Implementation**: The driver should support a `--convergence-check` mode
that automatically runs both migrations and reports a categorized diff
(code vs comments, algorithmic vs cosmetic).

### 3b. Divergence Patterns by Target

The dual-origin convergence test (§3) reveals different divergence counts
depending on the target mode.  Two target families are supported:

| Target | Type Declarations | Literal Mode | Intrinsic Mode |
|--------|-------------------|--------------|----------------|
| `kind16` | `REAL(KIND=16)` | `_16` suffix | `REAL(x, KIND=16)` |
| `multifloats` | `TYPE(float64x2)` | `float64x2(...)` constructor | Constructor wrapping |

**Why multifloats has fewer divergences than KIND-based targets:**

The multifloats target eliminates certain S/D asymmetries that KIND targets
preserve.  Two mechanisms drive this:

1. **Named-constant replacement.**  The `known_constants` map in the target
   YAML replaces local `PARAMETER` constants (`ONE`, `ZERO`, `HALF`, etc.)
   with module-provided values (`MF_ONE`, `MF_ZERO`, `MF_HALF`).  Many
   LAPACK D/Z routines define these constants locally while their S/C
   counterparts do not.  In KIND mode the D-half's `PARAMETER (ONE=1.0E0_16)`
   declaration survives into the output, creating a divergence.  In
   multifloats mode both halves use `MF_ONE` from `USE MULTIFLOATS`, so
   the local declaration is removed and the pair converges.

   *Example:* `dgehd2.f` defines `DOUBLE PRECISION ONE; PARAMETER(ONE=1.0D+0)`.
   `sgehd2.f` has no such constant.
   - KIND=16: D-half output has `REAL(KIND=16) ONE; PARAMETER(ONE=1.0E0_16)` → diverges (+2).
   - Multifloats: `ONE` is replaced by `MF_ONE` in both halves → converges.

2. **Constructor wrapping normalizes type conversions.**  In KIND mode,
   `REAL(x, KIND=16)` preserves the original intrinsic call structure, so
   differences in whether the S-half or D-half uses an explicit `REAL()`
   show up as divergences.  In multifloats mode, both `REAL(x)` and bare
   usage get wrapped uniformly in `float64x2(x)`, collapsing the
   distinction.

**Detailed divergence counts** are maintained in
`src-multifloats/DIVERGENCE.md`.

### 4. Compilation Test
Compile the migrated library with a Fortran compiler that supports `KIND=16`
(e.g., GFortran with `-freal-8-real-16`) and verify it compiles without errors.

### 5. Symbol Completeness
After migration, extract symbols from the migrated library and verify every
`D*` symbol has a corresponding `Q*` symbol.

### 6. Numerical Test
For a subset of routines, run the LAPACK test suite with the migrated library
and verify results match expected precision.
