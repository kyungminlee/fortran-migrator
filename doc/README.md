# fortran-migrator

**fortran-migrator** is an automated type-migration tool for Fortran and C codebases, built specifically for high-performance computing libraries such as BLAS, LAPACK, BLACS, PBLAS, and ScaLAPACK.

It parses existing fixed- and free-form Fortran source code to safely convert `REAL` and `COMPLEX` types from one kind to another, automatically updating the corresponding function names to match the new precision prefixes.

## Key Features
* **Automated Type Conversion:** Upgrades standard precision types to extended precisions (`KIND=10`, `KIND=16`) or multiword floating-point (`float64x2` via the multifloats target).
* **Smart Renaming:** Automatically updates function prefixes (e.g., converting `dgemm` to `qgemm` for KIND=16, or `ddgemm` for multifloats).
* **Format Agnostic:** Supports both fixed-form and free-form Fortran source code.
* **C Migration:** Template-based cloning for C libraries (BLACS, PBLAS) with mechanical type substitution.
* **Convergence Testing:** Dual-origin verification confirms migration correctness by comparing S-half and D-half migrated output.
* **YAML Recipes:** Declarative library descriptions drive the migration pipeline.

## Target Modes

### KIND-based targets

| Data Type | Target Extended Type | New Prefix | Example Conversion |
| :--- | :--- | :--- | :--- |
| `REAL` | `REAL(KIND=10)` * | `E` | `dgemm` → `egemm` |
| `REAL` | `REAL(KIND=16)` | `Q` | `dgemm` → `qgemm` |
| `COMPLEX` | `COMPLEX(KIND=10)` * | `Y` | `zgemm` → `ygemm` |
| `COMPLEX` | `COMPLEX(KIND=16)` | `X` | `zgemm` → `xgemm` |

*\* Note: `KIND=10` is specifically targeted for supported x86 architectures.*

### Multifloats target

Uses `float64x2` (double-double) arithmetic via an external module:

| Data Type | Target Type | New Prefix | Example Conversion |
| :--- | :--- | :--- | :--- |
| `REAL` | `TYPE(float64x2)` | `M` | `dgemm` → `mgemm` |
| `COMPLEX` | `TYPE(complex64x2)` | `W` | `zgemm` → `wgemm` |

## Implementation Details

Under the hood, **fortran-migrator** uses a **hybrid approach**: a compiler-based parser (Flang or GFortran) extracts structural facts from source files, and the Python engine applies regex-based transformations guided by those facts. This ensures syntactically aware conversions while preserving all formatting, comments, and preprocessor directives.

See [DEVELOPER.md](DEVELOPER.md) for the full developer guide.
