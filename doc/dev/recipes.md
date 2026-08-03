# Migration Recipes

A **recipe** is a YAML file that describes the structure of a numerical library and provides the parameters for its migration. The `migrator` uses these recipes to identify source files, discover symbols, and apply the correct transformations.

## Example: `codegen/recipes/blas.yaml`

```yaml
library: blas
language: fortran
source_dir: extern/lapack-3.12.1/BLAS/SRC
extensions: [.f]
```

Symbols are discovered by scanning `source_dir` for `SUBROUTINE`/`FUNCTION` definitions; the target prefix scheme comes from the target definition (`codegen/targets/*.yaml`), not the recipe.

## Recipe Fields

### Core Fields

#### `library` (string, required)
The short name of the library (e.g., `blas`, `lapack`, `scalapack`). Used for generating output filenames and library archives.

#### `language` (string, required)
The primary language of the source files:
*   `fortran`: Uses the hybrid parser-guided / regex migrator.
*   `c`: Uses the "clone-and-substitute" template migrator.

#### `source_dir` (string, required)
The relative path from the project root to the directory containing the source files.

#### `extensions` (list of strings, optional)
File extensions to include. Defaults to `[.f, .f90]`. Use `[.c, .h]` for C libraries. Note: `.F90` (uppercase) triggers preprocessing.

### File Handling

#### `skip_files` (list of strings, optional)
Source stems to skip during migration (case-insensitive). Useful for mixed-precision routines without a direct target equivalent (e.g., `DSDOT`).

#### `skip_files_manifest` (string, optional)
Recipe-relative file holding the same stems, one per line (`#` comments
and blank lines ignored). Unioned with `skip_files`. Use it for
generated bulk lists, so the recipe stays hand-readable and the list
lives next to its generator — `xblas.yaml` names
`xblas/skip_mixed_precision.txt`, the 402 mixed-input-precision stems
written by `recipes/xblas/gen_skip_list.py`.

#### `copy_files` (list of strings, optional)
Source stems to copy unchanged. Used for multi-precision utility routines that don't need transformation.

### Dependencies

#### `depends` (list of strings, optional)
Paths to dependency recipe files (resolved relative to the recipe's directory). Symbols from dependency libraries are included in the rename map. Example:
```yaml
depends:
  - blas.yaml
  - blacs.yaml
```

#### `extra_symbol_dirs` (list of strings, optional)
Extra directories to scan for symbols (resolved relative to project root). Used when symbols needed for renaming are defined outside the main `source_dir`:
```yaml
extra_symbol_dirs:
  - extern/lapack-3.12.1/INSTALL
```

### Precision Family Control

#### `prefer_source` (list of strings, optional)
Source stems whose S/C-migrated output should be kept as canonical instead of the default D/Z-first preference. Used when the S/C variant is more correct.

#### `local_renames` (dict, optional)
Maps S/C-half local identifiers to D/Z-half counterparts for convergence normalization:
```yaml
local_renames:
  CR: ZR
  SX: DX
```

### C-Specific Fields

#### `extra_c_dirs` (list of strings, optional)
Additional C source directories to migrate (flat-copied to output_dir). Resolved relative to project root.

#### `c_return_types` (list of strings, optional)
Regex fragments for additional C return types to recognize when scanning for function definitions:
```yaml
c_return_types:
  - 'PBTYP_T\s*\*'
```

#### `c_type_aliases` (list of dicts, optional)
Library-specific C typedef renames. Each entry maps a list of source identifiers to a target using template substitution:
```yaml
c_type_aliases:
  - from: [cmplx, cmplx8]
    to: 'cmplx{RPU}'
  - from: [cmplx16]
    to: 'cmplx{CPU}'
```
Template variables: `{REAL_TYPE}`, `{COMPLEX_TYPE}`, `{C_REAL_TYPE}`, `{RP}`, `{CP}`, `{RPU}`, `{CPU}`.

#### `header_patches` (list of dicts, optional)
Insert content into migrated headers. Each entry names the header
(`file`, relative to the migrated output directory), the text to insert
(`insert`, template-substituted), and one insertion point:
```yaml
header_patches:
  - file: TOOLS/PBblacs.h
    after: '#define MPI_DOUBLE_COMPLEX ...'
    insert: '#define MPI_QREAL ...'
    when: kind    # optional: only apply for this target family
```

Insertion points, checked in this order — `at_bof: true` (prepend at
file start), `at_eof: true` (append at file end), `before: <line>`
(insert on the line before the literal anchor), `after: <line>` (on the
line after it). An entry whose file is missing, whose anchor is absent,
or whose text is already present is skipped, so patching is idempotent.

`when` gates the entry on the target: `kind` (any KIND target),
`multifloats`, or an exact target name; absent means always applied.

### Overrides

#### `overrides` (dict, optional)
Target-gated verbatim file overrides. Hand-written files replace migrated output for a specific target:
```yaml
overrides:
  multifloats:
    src_dir: blacs/mfc_overrides
    files:
      - blacs_pinfo_.c
```

#### `patches` (list of strings, optional)
Unified-diff patch files (resolved relative to `codegen/recipes/<lib>/patches/`)
applied to the staged source tree before migration. The patched body
goes through the normal migration pipeline — so a single upstream-shape
patch produces correctly-renamed/promoted output for every target.
Used to carry bug fixes for upstream sources without editing
`extern/`.

```yaml
patches:
  - pdlanhs.f.patch
  - pzlanhs.f.patch
```

When a patch fixes a bug present only in the D/Z half (the migrator's
default canonical), pin the corresponding stem(s) in `prefer_source`
so the patched body wins convergence against the unpatched C/S
sibling — the canonical-rank picker doesn't recognize ScaLAPACK's
`pd*`/`pz*` two-letter prefix shape and would otherwise sort
`pclanhs.f` ahead of `pzlanhs.f`. See `doc/upstream-bugs/` for
the canonical examples.

#### `extra_renames` (dict, optional)
Force-rename entries appended to the classifier's rename map after
family discovery. Used for precision-prefixed orphan symbols whose
S/C sibling does not exist upstream (so the prefix classifier cannot
pair them into a precision family). Each entry maps an upstream
upper-cased identifier to a target template that may reference
`{RP}`/`{CP}`/`{RPU}`/`{CPU}` via target template_vars.

```yaml
extra_renames:
  PDLAIECTB: P{RPU}LAIECTB
  PDLAIECTL: P{RPU}LAIECTL
```

The canonical example is ScaLAPACK's `pdlaiectb_`/`pdlaiectl_`
helpers, which exist only in double precision because the bit-shift
sign-bit hack they rely on is 64-bit-only. Without the rename, the
migrated clones would still export `pdlaiectb_`/`pdlaiectl_`,
clashing with the standard archive's identically-named symbols.

#### `ep_renames` (dict, optional)
Rename calls into the `ep_`-privatized support engine, resolved against
the privatization manifest instead of spelled out by hand:
```yaml
ep_renames:
  manifest: privatize_ep_symbols.txt
  symbols:
    - BLACS_GRIDINIT
    - DESCINIT
```
Each name becomes an `extra_renames` entry `NAME → EP_NAME` — the
source-level form of the `name → ep_name` rewrite `privatize.py` applies
to the recipes that carry `privatize_symbols:`. A name the manifest does
not carry is a load error.

For a recipe that calls into the privatized engine without being
privatized itself (MUMPS is the only one), this is what ties the two
together: dropping a symbol from the manifest then fails the recipe
load, rather than surfacing as an undefined `ep_<name>_` at final
executable link, far from its cause.
