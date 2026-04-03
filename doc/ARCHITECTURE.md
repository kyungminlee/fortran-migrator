# Architecture

**fortran-migrator** is an automated source-to-source translation tool for numerical libraries. Its design handles the unique challenges of Fortran (fixed-form column width, continuation lines) and high-performance computing (maintaining identical source layout, preserving preprocessor macros).

## Hybrid Migration Strategy

The tool uses a **hybrid approach** that combines the accuracy of a compiler-based parse tree with the formatting preservation of source-level rewriting.

1.  **Parse with Flang**: The tool invokes `flang-new` (the LLVM Fortran compiler) as a subprocess to generate a textual dump of the source's parse tree.
2.  **Extract Facts**: The `pyengine/flang_parser.py` scans the parse tree dump to extract key information (facts) about the source file:
    *   **Type Declarations**: Location and nature of `REAL`, `DOUBLE PRECISION`, etc.
    *   **Routine Definitions**: `SUBROUTINE` and `FUNCTION` names.
    *   **Call Sites**: Locations where other routines are called.
    *   **Literals**: Floating-point constants with D or E exponents.
    *   **Intrinsics**: Calls to type-specific intrinsic functions.
3.  **Apply Transformations**: The `pyengine/fortran_migrator.py` uses these facts as an "oracle" to guide line-by-line regex replacements. This ensures that only relevant code is modified while comments, whitespace, and preprocessor directives are left exactly as-is.
4.  **Reformat**: For fixed-form Fortran, the tool automatically reformats lines that exceed the 72-character limit after transformation, ensuring they remain syntactically valid.

## Core Components

### Python Engine (`pyengine/`)
The primary driver of the tool. It orchestrates the migration pipeline:
*   **`pipeline.py`**: Manages the multi-file migration process.
*   **`symbol_scanner.py`**: Discovers the API of the source library to build the rename map.
*   **`fortran_migrator.py`**: Implements the Fortran rewriting logic.
*   **`c_migrator.py`**: Implements the C template-based cloning (used for BLACS).
*   **`flang_parser.py`**: Interfaces with the Flang compiler.

### Symbol Database
The migrator uses a **Symbol Database** to drive naming decisions. By identifying which routines follow the standard precision prefix convention (S/D/C/Z), it can safely rename `dgemm` to `qgemm` while leaving type-independent routines like `LSAME` or `XERBLA` untouched.

### C++ Migrator (`src/`)
A high-performance C++ implementation is planned (currently in the stub phase). The goal is to eventually replace the Python engine's rewriting logic with a direct C++ library that links against `libflangParser`. This will provide:
*   **Direct AST access**: Eliminating the need to parse textual parse tree dumps.
*   **Byte-level rewriting**: Using exact source offsets for even more precise transformations.
*   **Standalone binary**: A single, fast executable with no Python dependencies.

## C Migration (BLACS Style)

For C-based libraries like BLACS, the tool uses a **template-based cloning** strategy. Since these libraries are often written with nearly identical files for each precision, the migrator:
1.  Identifies a "template" file (e.g., `dgebr2d_.c`).
2.  Performs mechanical text substitution on C types (e.g., `double` → `QREAL`).
3.  Updates MPI datatype constants (e.g., `MPI_DOUBLE` → `MPI_QREAL`).
4.  Generates a new file with the correct prefix (e.g., `qgebr2d_.c`).

This process is handled by `pyengine/c_migrator.py`.
