# Migration Recipes

A **recipe** is a YAML file that describes the structure of a numerical library and provides the parameters for its migration. The `pyengine` uses these recipes to identify source files, discover symbols, and apply the correct transformations.

## Example: `recipes/blas.yaml`

```yaml
library: blas
language: fortran
source_dir: external/lapack-3.12.1/BLAS/SRC
extensions: [.f]

symbols:
  method: scan_source

prefix:
  style: direct
```

## Recipe Fields

### `library` (string, required)
The short name of the library (e.g., `blas`, `lapack`, `scalapack`). This is used for generating output filenames and library archives.

### `language` (string, required)
The primary language of the source files. Supported values:
*   `fortran`: Uses the hybrid Flang-guided / regex migrator.
*   `c`: Uses the "clone-and-substitute" template migrator (used for BLACS).

### `source_dir` (string, required)
The relative path from the project root to the directory containing the source files.

### `extensions` (list of strings, optional)
The list of file extensions to include in the migration. Defaults to `[.f, .f90]`.

### `symbols` (object, optional)
Configuration for how routine names are discovered for the rename map.
*   `method`:
    *   `scan_source` (default): Scans the `source_dir` for `SUBROUTINE` and `FUNCTION` definitions.
    *   `nm_library`: Uses the `nm` tool to extract symbols from a pre-built static or shared library.
*   `library_path`: The path to the pre-built library file (required if `method` is `nm_library`).

### `prefix` (object, optional)
Configuration for the naming convention style.
*   `style`:
    *   `direct` (default): Single-character prefix (e.g., `DGEMM` → `QGEMM`).
    *   `scalapack`: `P` + precision character (e.g., `PDGESV` → `PQGESV`).

### `skip_files` (list of strings, optional)
A list of filenames (without extension) to skip during migration. This is useful for:
*   Mixed-precision routines that don't have a direct target equivalent (e.g., `DSDOT`).
*   Iterative-refinement routines that require manual porting.

### `copy_all_originals` (boolean, optional)
*For C migration only.* If `true`, all original source files are copied to the output directory before the migrated "clones" are added. This is useful for libraries where only a subset of files are precision-specific.

### `patches` (list of strings, optional)
*Future work.* A list of patch files to apply to the migrated source to handle complex manual adjustments.

## Symbol Discovery Strategies

### `scan_source`
The tool scans the `source_dir` for Fortran definitions. It identifies which files define which routines and uses this to build the prefix map.
*   **Pros**: No pre-built library required.
*   **Cons**: May miss symbols hidden behind preprocessor directives or non-standard syntax.

### `nm_library`
The tool runs `nm --defined-only` on the provided library file.
*   **Pros**: Most accurate way to find the actual public API of the library.
*   **Cons**: Requires the library to be compiled in its original precision first.
