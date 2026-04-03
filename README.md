# fortran-migrator

**fortran-migrator** is an automated type-migration tool for Fortran codebases, built specifically for high-performance computing libraries such as BLAS, LAPACK, BLACS, and ScaLAPACK. 

It parses existing fixed- and free-form Fortran source code to safely convert `REAL` and `COMPLEX` types from one kind to another, automatically updating the corresponding function names to match the new precision prefixes.

## Key Features
* **Automated Type Conversion:** Seamlessly upgrades standard precision types to extended precisions (`KIND=10` and `KIND=16`).
* **Smart Renaming:** Automatically updates function prefixes (e.g., converting `dgemm` to `qgemm`, or `zgemm` to `xgemm`).
* **Format Agnostic:** Supports both fixed-form and free-form Fortran source code.
* **Macro Preservation:** Safely retains all preprocessor macros during the parse and rewrite cycle.

## Type Mapping & Naming Conventions

In standard BLAS/LAPACK, function prefixes denote the data type (`s` for `REAL`, `d` for `DOUBLE PRECISION`, `c` for `COMPLEX`, and `z` for `DOUBLE COMPLEX`). 

**fortran-migrator** extends this convention to support 80-bit and 128-bit floating-point types:

| Data Type | Target Extended Type | New Prefix | Example Conversion |
| :--- | :--- | :--- | :--- |
| `REAL` | `REAL(KIND=10)` * | `e` | `dgemm` → `egemm` |
| `REAL` | `REAL(KIND=16)` | `q` | `dgemm` → `qgemm` |
| `COMPLEX` | `COMPLEX(KIND=10)` * | `y` | `zgemm` → `ygemm` |
| `COMPLEX` | `COMPLEX(KIND=16)` | `x` | `zgemm` → `xgemm` |

*\* Note: `KIND=10` is specifically targeted for supported x86 architectures.*

## Implementation Details

Under the hood, **fortran-migrator** leverages the **LLVM/Flang** compiler infrastructure. By generating and inspecting the parse tree of the Fortran source files, the tool guarantees syntactically aware conversions rather than relying on simple text replacement. This allows it to accurately identify type declarations and function calls while leaving preprocessor directives completely intact.
It also takes symbols from a built library (static or shared) for smart renaming.
