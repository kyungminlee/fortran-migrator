# Architecture

**fortran-migrator** is an automated source-to-source translation tool for numerical libraries. Its design handles the unique challenges of Fortran (fixed-form column width, continuation lines) and high-performance computing (maintaining identical source layout, preserving preprocessor macros).

## Hybrid Migration Strategy

The tool uses a **hybrid approach** that combines the accuracy of a compiler-based parse tree with the formatting preservation of source-level rewriting.

1.  **Parse with a Fortran compiler**: The tool invokes either `flang-new` (LLVM) or `gfortran` as a subprocess to generate a textual dump of the source's parse tree. The `--parser` flag selects the backend; without it, the tool falls back to regex-only scanning.
2.  **Extract Facts**: The parser module (`pyengine/flang_parser.py` or `pyengine/gfortran_parser.py`) scans the parse tree dump to extract key information:
    *   **Type Declarations**: Location and nature of `REAL`, `DOUBLE PRECISION`, etc.
    *   **Routine Definitions**: `SUBROUTINE` and `FUNCTION` names.
    *   **Call Sites**: Locations where other routines are called.
    *   **Literals**: Floating-point constants with D or E exponents.
    *   **Intrinsics**: Calls to type-specific intrinsic functions.
3.  **Apply Transformations**: The `pyengine/fortran_migrator.py` uses these facts as an "oracle" to guide line-by-line regex replacements. This ensures that only relevant code is modified while comments, whitespace, and preprocessor directives are left exactly as-is.
4.  **Reformat**: For fixed-form Fortran, the tool automatically reformats lines that exceed the 72-character limit after transformation, ensuring they remain syntactically valid.

## Core Components

### Python Engine (`src-multifloats/pyengine/`)
The primary driver of the tool. It orchestrates the migration pipeline:
*   **`pipeline.py`**: Manages the multi-file migration process.
*   **`target_mode.py`**: Defines `TargetMode` — the abstraction for different target precisions (KIND=10/16, multifloats).
*   **`config.py`**: Loads and validates YAML recipe files.
*   **`symbol_scanner.py`**: Discovers the API of the source library by scanning source files or compiled libraries.
*   **`prefix_classifier.py`**: Classifies routine names by precision prefix and builds the rename map.
*   **`fortran_migrator.py`**: Implements the Fortran rewriting logic.
*   **`c_migrator.py`**: Implements the C template-based cloning (used for BLACS, PBLAS).
*   **`flang_parser.py`**: Interfaces with the Flang compiler for parse tree extraction.
*   **`gfortran_parser.py`**: Interfaces with GFortran as an alternative parser backend.

### Symbol Database
The migrator uses a **prefix classifier** to drive naming decisions. By scanning source files (or using `nm` on compiled libraries) to identify which routines follow the standard precision prefix convention (S/D/C/Z), it can safely rename `dgemm` to `qgemm` while leaving type-independent routines like `LSAME` or `XERBLA` untouched.

## C Migration (BLACS / PBLAS Style)

For C-based libraries, the tool uses a **template-based cloning** strategy. Since these libraries are written with nearly identical files for each precision, the migrator:
1.  Identifies a "template" file (e.g., `dgebr2d_.c`).
2.  Performs mechanical text substitution on C types (e.g., `double` → `__float128` or `float64x2_t`).
3.  Updates MPI datatype constants (e.g., `MPI_DOUBLE` → `MPI_FLOAT128`).
4.  Generates a new file with the correct prefix (e.g., `qgebr2d_.c`).

This process is handled by `pyengine/c_migrator.py`.
