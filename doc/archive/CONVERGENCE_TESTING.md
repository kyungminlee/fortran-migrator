# Convergence Testing

**Convergence Testing** is a powerful verification technique used to ensure the accuracy and completeness of the migration process. It leverages the dual-origin nature of numerical libraries like BLAS and LAPACK.

## The Principle

Almost every BLAS/LAPACK routine exists in at least two variants:
*   **Single Precision**: `S` (real) or `C` (complex)
*   **Double Precision**: `D` (real) or `Z` (complex)

When both the single-precision and double-precision variants of a routine are migrated to the same target precision (e.g., `KIND=16`), the resulting source code should be **identical** (modulo algorithmic differences or comments).

```
sgemm.f ──migrate(KIND=16)──→ qgemm_from_s.f ─┐
                                                ├── diff
dgemm.f ──migrate(KIND=16)──→ qgemm_from_d.f ─┘
```

## Why It Works

The migrator is designed to normalize all precision-dependent constructs:
*   **Type Declarations**: `REAL` and `DOUBLE PRECISION` both become `REAL(KIND=16)`.
*   **Literals**: `1.0E+0` and `1.0D+0` both become `1.0E+0_16`.
*   **Intrinsics**: `REAL(x)` and `DBLE(x)` both become `REAL(x, KIND=16)`.
*   **Routine Names**: `SGEMM` and `DGEMM` both become `QGEMM`.

If the two migrated outputs are identical, it confirms that the migrator has correctly handled all precision-specific variants on both paths.

## What Mismatches Reveal

If the migrated files differ, the `diff` provides valuable diagnostic information:

### 1. Algorithmic Differences
Some routines are implemented differently for single and double precision (e.g., different iteration counts, different convergence thresholds). These are **intentional** and the diff highlights them for manual review to determine which approach is better for quad precision.

### 2. Hardcoded Constants
A literal `1.0E-6` in a single-precision routine vs. `1.0D-12` in the double-precision version indicates a precision-dependent tolerance. The migrator cannot automatically determine the "correct" tolerance for quad precision, so these values must be manually adjusted in the target code.

### 3. Missing Conversions
If the `diff` shows a `DBLE()` call that wasn't converted, or a routine name that didn't get its prefix updated, it reveals a bug or an unhandled case in the migrator's logic.

## Performing the Test

To perform a convergence check manually:
1.  Migrate the single-precision source to the target KIND.
2.  Migrate the double-precision source to the same target KIND.
3.  Use a diff tool (e.g., `diff -u`) to compare the results.

Future versions of `fortran-migrator` will include an automated `--convergence-check` flag that performs this process for all available routine pairs and generates a summary report.
