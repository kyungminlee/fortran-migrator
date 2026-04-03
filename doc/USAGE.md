# Usage Guide

**fortran-migrator** provides a multi-step pipeline to migrate numerical Fortran/C libraries to extended precision (`KIND=10` or `KIND=16`).

The tool is primarily used through the Python engine (`pyengine`).

## Prerequisite: uv

This project uses `uv` for Python dependency management. Ensure you have it installed:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

## Basic Workflow

The most common way to use the tool is the `run` command, which executes the full pipeline:
1.  **Migrate**: Rewrite source files to the target precision.
2.  **Verify**: Check the output for residual types, literals, or column-width overflows.
3.  **Compile**: Verify the syntax of the migrated code using a Fortran compiler.
4.  **Build**: Create static libraries from the migrated source.

```bash
uv run python -m pyengine run recipes/blas.yaml work/ --kind 16
```

## CLI Commands

### `migrate`
Only performs the source-to-source rewriting step.

```bash
uv run python -m pyengine migrate recipes/lapack.yaml output/ --kind 16
```
*   `recipe`: Path to the library's YAML recipe.
*   `output_dir`: Where to write the migrated files.
*   `--kind`: Target precision (10 for 80-bit, 16 for 128-bit).
*   `--dry-run`: Show what would be changed without writing files.

### `verify`
Performs a set of heuristic checks on the migrated source files to identify potential issues.

```bash
uv run python -m pyengine verify output/
```
Checks for:
*   Residual precision types (e.g., `DOUBLE PRECISION` that wasn't converted).
*   Residual `D`-exponent literals (e.g., `1.0D0`).
*   Column-width overflows in fixed-form code (lines exceeding 72 characters).

### `compile`
Attempts to compile every migrated file to verify syntactic correctness.

```bash
uv run python -m pyengine compile output/ --fc gfortran
```
*   `--fc`: Path to the Fortran compiler (defaults to `gfortran`).

### `build`
Compiles and archives the migrated files into static libraries. It produces two libraries:
1.  `lib<prefix><library>.a`: Contains precision-specific routines (e.g., `libqblas.a`).
2.  `lib<library>_common.a`: Contains precision-independent routines.

```bash
uv run python -m pyengine build recipes/blas.yaml output/ --kind 16 --fc gfortran
```

### `run`
Runs all the above steps in sequence.

```bash
uv run python -m pyengine run recipes/blas.yaml work/ --kind 16
```

## Target Precisions

The tool supports two target precision "KIND" values:

| KIND | Description | Precision Prefix (Real/Complex) |
| :--- | :--- | :--- |
| **10** | 80-bit extended precision (x86) | `E` / `Y` |
| **16** | 128-bit quad precision | `Q` / `X` |

Example: If migrating BLAS to `KIND=16`, `dgemm` becomes `qgemm` and `zgemm` becomes `xgemm`.
