# Review & Revised Plan: `multifloats` Support

## 1. Assessment of the Original Plan

### Strengths

- **Option D (Semantic Oracle + Source Rewriting) is the correct architecture.**
  The existing codebase already implements this pattern: `flang_parser.py` extracts
  `ParseTreeFacts` and the regex engine in `fortran_migrator.py` applies line-level
  transformations. The plan extends this rather than inventing a new paradigm.

- **The literal conversion strategy is sound.** Distinguishing
  `float64x2('1.0D+0')` (string constructor, lossless) from `float64x2(1.0D0)`
  (double constructor, for exact values) avoids a subtle class of bugs where an
  intermediate double representation silently loses precision.

- **The PARAMETER/DATA strategy is realistic.** The two-tiered approach for
  PARAMETER (module-imported constants for common values, runtime-initialized
  variables for uncommon ones) correctly handles the Fortran standard restriction
  that derived-type constructors are not constant expressions.

- **The transformation difference table (Section "Key Differences") is precise.**
  It identifies the six dimensions where multifloats diverges from KIND-based
  migration — a critical scoping artifact.

- **Phase 0 (prototype with `dgemv`) is prudent.** Manual migration of a real
  routine before writing automation catches unknown-unknowns early.

- **Each phase maps cleanly to existing code.** Phase 2 maps to
  `replace_type_decls()`, Phase 4 to `replace_literals()`, Phase 7 to
  `replace_intrinsic_calls()`, etc.

### Weaknesses

- The plan assumes `kind` remains an integer throughout the codebase, but
  multifloats is not a Fortran KIND — it is a fundamentally different target mode.

- No mention of updating the convergence/divergence checking machinery, which is
  tightly coupled to KIND-based output syntax.

- No mention of recipe configuration changes (`RecipeConfig` fields for module
  name, type names, constructor names, known-constants map).

- No mention of `la_constants` module rewriting, which is significant for LAPACK.

- Heavy focus on fixed-form BLAS; free-form LAPACK `.f90` files that use
  `integer, parameter :: wp = kind(1.d0)` are not addressed.

---

## 2. Assessment of Gemini's Critique

### Valid Points

- **Phase 8 is scheduled too late.** The semantic oracle is needed to correctly
  implement type declarations (Phase 2) and literal replacement (Phase 4). However,
  moving it all the way to Phase 1 creates a chicken-and-egg problem (need compiled
  multifloats module for sema). Phase 1.5 is the right position.

- **SAVE semantics in DATA conversion.** In Fortran, `DATA` implies `SAVE` — the
  variable retains its value across calls. Converting `DATA X / 1.0D0 /` to a plain
  executable assignment `X = float64x2('1.0')` means X is re-initialized on every
  call. The migrator must add an explicit `SAVE` attribute. The original plan missed
  this.

- **COMMON and EQUIVALENCE blocks.** Derived types in COMMON require `SEQUENCE`.
  For BLAS/LAPACK specifically, COMMON blocks with FP members do not exist and
  EQUIVALENCE is limited to 2 files (`dlaln2.f`, `slaln2.f`), so this is a
  "flag and warn" issue rather than a blocker.

- **Expanding Phase 0 prototype scope.** `dgemv` lacks DATA statements. Adding
  `drotm.f` (which has `DATA ZERO,TWO/0.D0,2.D0/`) exercises the hard
  transformations early.

### Overstated Points

- **"Module Mocking" as a separate phase.** The original plan already addresses this
  in Phase 0 ("Create a multifloats stub module") and Phase 8 (fallback to
  no-sema + heuristics). A stub `.f90` module is a one-time task inside Phase 0.

- **"Column-Aware Rewriter" referencing `rewriter.cpp`.** No such file exists.
  The actual logic is `reformat_fixed_line()` in `fortran_migrator.py:556`. The
  concern about column overflow is valid but the solution is already partially in
  place — it needs enhancement, not a new component.

- **Moving the semantic oracle to Phase 1 is too aggressive.** It creates a
  dependency on having a compiled multifloats `.mod` file before any infrastructure
  exists. Phase 1.5 is more pragmatic.

- **PARAMETER-in-constant-expression context** is theoretically valid but
  practically minor. In BLAS, PARAMETER constants (ZERO, ONE) are used exclusively
  in executable expressions, not as array bounds. Still worth a diagnostic check.

  One notable exception: `dlaruv.f` contains
  `PARAMETER ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )` where `R` is defined in
  terms of another PARAMETER — these chained dependencies need careful handling.

---

## 3. Gaps in Both Documents

### 3.1 The `kind` Parameter is an `int` Everywhere

In `__main__.py`, `--kind` is `type=int, choices=[10, 16]`. In `pipeline.py`,
`run_migration()` takes `target_kind: int`. In `fortran_migrator.py`, every function
signature uses `kind: int`. Multifloats is not an integer kind — it requires either
a `str | int` discriminated union or a new `TargetMode` enum/dataclass. This is the
single most pervasive change to the codebase.

### 3.2 Convergence/Divergence Pipeline

`pipeline.py` has elaborate convergence checking: `_canonicalize_for_compare()`,
`_light_normalize()`, `_strip_real_cmplx_casts()`, `_filter_precision_drift()`.
These assume KIND-based output syntax. For multifloats, they must handle
`float64x2(...)` constructors, `TYPE(float64x2)` declarations, and
`USE multifloats` statements.

### 3.3 Recipe Configuration

`RecipeConfig` needs new fields for multifloats: module name, real/complex type
names, constructor names, and the known-constants map (ZERO -> MF_ZERO, etc.).
These must be configurable per recipe.

### 3.4 Prefix Map and Rename Map

`PREFIX_MAP` in `prefix_classifier.py` only knows `{10, 16}` as kind targets.
This needs an entry for multifloats or a bypass mechanism. Routine renaming (D ->
Q/E/X) is orthogonal to the type system change and should still work, but the
`build_rename_map()` code path must accept the new target mode.

### 3.5 Free-Form `.f90` Files

LAPACK's newer files use `integer, parameter :: wp = kind(1.d0)`. The existing
migrator handles these via `_replace_kind_parameter()`. For multifloats, `wp` cannot
be an integer KIND parameter — it must be removed entirely and replaced with direct
`TYPE(float64x2)` usage. This is a substantial transformation not discussed in
either document.

### 3.6 `la_constants` Module

`fortran_migrator.py:680-768` implements `rewrite_la_constants_use()` which renames
`LA_CONSTANTS` to `LA_CONSTANTS_EP` and remaps precision-tagged constants
(SZERO -> QZERO, etc.). For multifloats, the constants module needs to provide
`float64x2` constants instead. This interaction is undiscussed.

### 3.7 The `gfortran_parser.py` Backend

The codebase supports both `flang` and `gfortran` parser backends. The plan only
mentions flang. If gfortran's parse tree dump cannot provide the needed semantic
information, multifloats may be flang-only.

### 3.8 Build System

CMake generation in `__main__.py` (`_generate_cmake()`) assumes standard Fortran
compilation. Multifloats requires linking against the multifloats module — the
CMake template needs a multifloats-aware variant with `find_package` or
`target_link_libraries`.

---

## 4. Technical Risks

### 4.1 Multi-Line Statement Processing

The regex engine processes one line at a time (`migrate_fixed_form()` at line 597).
DATA and PARAMETER statements can span continuation lines:
```fortran
      DATA (A(I), I=1,10)
     + / 1.0D0, 2.0D0, 3.0D0,
     + 4.0D0, 5.0D0, 6.0D0,
     + 7.0D0, 8.0D0, 9.0D0, 10.0D0 /
```
Transforming these requires joining continuations, parsing the full statement, then
re-emitting. The existing `_dedup_intrinsic_stmts()` does multi-line joining for
INTRINSIC statements, but there is no equivalent for DATA. This is non-trivial.

### 4.2 Self-Interference in Literal Replacement

`replace_literals()` uses regex on code segments after splitting out string literals.
For multifloats, the replacement `float64x2('1.0D+0')` *introduces* a string literal
into the code. If the function is called in a second pass, the newly-introduced
strings must not be re-processed. The current string-splitting logic would protect
them, but this needs careful verification.

### 4.3 Fixed-Form Line Breaking Inside Constructors

`reformat_fixed_line()` does naive 66-character chunking, breaking at commas and
spaces. `float64x2('1.0D+0')` is 21 characters and must not be broken
mid-constructor. The break-point logic needs to avoid splitting inside
`float64x2(...)` calls.

### 4.4 Intrinsic Replacement Needs a Third Mode

`replace_intrinsic_calls()` has two modes: `needs_kind=False` (simple rename) and
`needs_kind=True` (add KIND argument). For multifloats, conversion intrinsics like
DBLE(x) become `float64x2(x)` — neither a simple rename nor a KIND addition. A
third replacement mode is needed: "replace with constructor call."

### 4.5 PARAMETER Chaining

`dlaruv.f` contains:
```fortran
      DOUBLE PRECISION ONE
      PARAMETER ( ONE = 1.0D0 )
      DOUBLE PRECISION R
      PARAMETER ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )
```
`R` depends on `ONE`. If `ONE` is imported from the multifloats module and `R` is
converted to a runtime variable, the initialization order matters. The migrator must
detect and correctly order such chains.

### 4.6 EQUIVALENCE Blocks

`dlaln2.f` and `slaln2.f` contain EQUIVALENCE statements mapping 2D arrays to 1D
vectors. If the equivalenced variables are floating-point and migrated to derived
types, EQUIVALENCE becomes invalid (Fortran standard prohibits EQUIVALENCE with
non-SEQUENCE derived types). These must be flagged for manual review.

---

## 5. Revised Implementation Plan

### Phase 0: Manual Prototype

**Goal:** End-to-end manual migration to validate the transformation catalog.

**Targets:**
- `dgemv.f` — clean routine, validates type decls + literals + intrinsics + USE
- `drotm.f` — has `DATA ZERO,TWO/0.D0,2.D0/`, validates DATA conversion + SAVE
- `drotmg.f` — has DATA with non-trivial values (GAM, GAMSQ, RGAMSQ)
- One PARAMETER-heavy BLAS file (e.g., `dgemm.f` with `PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)`)

**Tasks:**
- [ ] Hand-migrate each file to use `float64x2` / `multifloats`
- [ ] Catalog every transformation applied (type, line, before/after)
- [ ] Create a `multifloats` stub module (`.f90`) sufficient for compilation
- [ ] Verify the stub compiles and migrated files compile against it
- [ ] Identify any transformation categories not covered in this plan

### Phase 1: Target Mode Infrastructure

**Goal:** The migrator accepts `--target multifloats` and activates derived-type
transformation paths.

**Key design decision:** Introduce a `TargetMode` abstraction that encapsulates
the target kind. This replaces the raw `int` kind parameter throughout the codebase.

```python
@dataclass
class TargetMode:
    name: str            # 'kind10', 'kind16', 'multifloats'
    kind: int | None     # 10, 16, or None for multifloats
    real_type: str       # 'REAL(KIND=10)', 'TYPE(float64x2)'
    complex_type: str    # 'COMPLEX(KIND=16)', 'TYPE(complex128x2)'
    module_name: str | None          # None for KIND, 'multifloats' for MF
    real_constructor: str | None     # 'float64x2'
    complex_constructor: str | None  # 'complex128x2'
    known_constants: dict[str, str]  # {'ZERO': 'MF_ZERO', 'ONE': 'MF_ONE', ...}
```

**Files affected:**
- `pyengine/__main__.py` — new `--target` CLI option (or extend `--kind`)
- `pyengine/config.py` — new `RecipeConfig` fields for multifloats config
- `pyengine/pipeline.py` — propagate `TargetMode` instead of `int`
- `pyengine/fortran_migrator.py` — accept `TargetMode` in all transform functions
- `pyengine/prefix_classifier.py` — extend `PREFIX_MAP` for multifloats

**Tasks:**
- [ ] Design and implement `TargetMode` dataclass
- [ ] Update CLI to accept `--target multifloats`
- [ ] Update `RecipeConfig` with multifloats configuration fields
- [ ] Propagate `TargetMode` through the pipeline (replace `kind: int`)
- [ ] Gate every transformation function behind a mode check
- [ ] Ensure KIND=10/16 paths are completely unaffected (regression safety)

### Phase 1.5: Semantic Oracle Enrichment

**Goal:** Extract richer semantic facts from `flang-new` with sema enabled.

**Rationale for this position:** Phases 2-6 all benefit from semantic information.
However, this phase requires the multifloats stub module from Phase 0, so it cannot
come earlier. It must precede Phase 2 because type declaration rewriting needs to
know which variables are floating-point.

**Files affected:**
- `pyengine/flang_parser.py` — sema-enabled invocation and new fact extraction

**New `ParseTreeFacts` fields:**
- `parameter_names: list[tuple[str, str]]` — (name, value) for FP PARAMETERs
- `data_stmt_vars: list[tuple[str, str]]` — (name, value) for FP DATA entries
- `procedure_boundaries: list[tuple[str, int, int]]` — (name, start_line, end_line)
- `executable_boundary: dict[str, int]` — procedure_name -> first executable line

**Tasks:**
- [ ] Add `use_sema: bool` parameter to `run_flang_parse_tree()`
- [ ] When `use_sema=True`, use `-fdebug-dump-parse-tree` (without `-no-sema`)
- [ ] Provide the multifloats stub `.mod` file path via `-I` flag
- [ ] Extract PARAMETER context facts (which names, which values)
- [ ] Extract DATA statement facts (which variables, which values)
- [ ] Extract procedure boundaries and executable section start lines
- [ ] Handle sema failures with fallback to no-sema + heuristics
- [ ] Determine whether `gfortran_parser.py` can provide equivalent information

### Phase 2: Type Declaration Rewriting

**Goal:** `DOUBLE PRECISION` / `REAL(KIND=n)` -> `TYPE(float64x2)`,
`COMPLEX*16` / `COMPLEX(KIND=n)` -> `TYPE(complex128x2)`.

**Files affected:**
- `pyengine/fortran_migrator.py` — parallel function `replace_type_decls_multifloats()`

**Implementation notes:**
- The output syntax differs fundamentally: `TYPE(float64x2)` vs `REAL(KIND=16)`.
  This requires new regex patterns, not just different substitution strings.
- `TYPE(float64x2)` is 15 characters, comparable to `DOUBLE PRECISION` (16 chars),
  so column overflow is manageable for declarations.
- Must handle all input forms: `DOUBLE PRECISION`, `REAL*8`, `REAL(KIND=8)`,
  `REAL(8)`, `DOUBLE COMPLEX`, `COMPLEX*16`, `COMPLEX(KIND=8)`.
- Must handle function return types:
  `DOUBLE PRECISION FUNCTION FOO` -> `TYPE(float64x2) FUNCTION FOO`.
- Must preserve attributes: `INTENT`, `DIMENSION`, `ALLOCATABLE`, etc.
  Note: `TYPE(float64x2), INTENT(IN) :: X` is valid Fortran.

**Tasks:**
- [ ] New regex patterns for `TYPE(float64x2)` output
- [ ] Handle all input declaration forms
- [ ] Handle function return type declarations
- [ ] Preserve all attributes
- [ ] Handle `replace_standalone_real_complex()` equivalent for multifloats

### Phase 3: USE Statement Insertion

**Goal:** Every procedure using `float64x2` gets `USE multifloats`.

**Files affected:**
- `pyengine/fortran_migrator.py` — new function `insert_use_multifloats()`

**Implementation notes:**
- Must detect procedure boundaries: `SUBROUTINE`, `FUNCTION`, `PROGRAM`,
  `MODULE` (though MODULE is rare in BLAS/LAPACK fixed-form code).
- Insert position: immediately after the procedure header line
  (after `SUBROUTINE FOO(A, B, C)`), before `IMPLICIT NONE` and other declarations.
- In fixed-form: the USE statement must start in column 7 and fit within column 72.
  `      USE multifloats` is 22 characters — fits easily.
- The ONLY clause should list needed constants when PARAMETER lines are removed
  (coordinate with Phase 5):
  `USE multifloats, ONLY: float64x2, ZERO => MF_ZERO, ONE => MF_ONE`
- Avoid duplicate insertion if `USE multifloats` already present.

**Tasks:**
- [ ] Detect procedure boundaries (first pass or from `ParseTreeFacts`)
- [ ] Insert `USE multifloats` after procedure header
- [ ] Build ONLY clause dynamically based on needed imports
- [ ] Avoid duplicates
- [ ] Handle fixed-form column constraints

### Phase 4: Literal Replacement

**Goal:** `1.0D+0` -> `float64x2('1.0D+0')` or `float64x2(1.0D0)`.

**Files affected:**
- `pyengine/fortran_migrator.py` — new function `replace_literals_multifloats()`

**Implementation notes:**
- Constructor calls are significantly longer than KIND suffixes:
  `1.0D+0` (6 chars) -> `float64x2('1.0D+0')` (21 chars).
  This will trigger column-72 overflows in fixed-form code. The line-breaking
  engine must be enhanced to avoid splitting inside constructor calls.
- Classification of exact-in-double values: `0.0`, `1.0`, `2.0`, `-1.0`, `0.5`,
  powers of 2, etc. These can use the shorter `float64x2(1.0D0)` form.
- Must NOT replace literals inside PARAMETER statements (handled in Phase 5).
- Must NOT replace integer literals or character literals.
- Must NOT replace literals inside string literals (the current string-splitting
  logic handles this, but verify with the new constructor syntax which itself
  contains string literals).

**Tasks:**
- [ ] New literal replacement logic emitting constructor calls
- [ ] Classify literals as exact-in-double vs. requires-string-constructor
- [ ] Skip literals inside PARAMETER/DATA statements (use `ParseTreeFacts`)
- [ ] Enhance `reformat_fixed_line()` to avoid breaking inside constructor calls
- [ ] Verify no self-interference (constructor-introduced strings not re-processed)

### Phase 5: PARAMETER Statement Conversion

**Goal:** Remove `PARAMETER` declarations of FP constants; replace with module
imports or runtime initialization.

**Files affected:**
- `pyengine/fortran_migrator.py` — new function `convert_parameter_stmts()`
- `pyengine/flang_parser.py` — PARAMETER context facts (from Phase 1.5)

**Implementation notes:**
- Known constants (configurable per recipe via `known_constants` map):
  Remove PARAMETER line, add name to USE..ONLY import list, remove old type
  declaration.
- Unknown constants: Convert to `TYPE(float64x2) :: name` declaration +
  assignment in executable section. Must detect the executable section boundary.
- Chained PARAMETERs (e.g., `R = ONE / IPW2` in `dlaruv.f`): must detect
  dependency chains and order assignments correctly.
- PARAMETERs used in constant-expression contexts (array bounds, other PARAMETER
  initializers): must be flagged — these cannot be converted to runtime variables.

**Tasks:**
- [ ] Parse PARAMETER statements (multi-line, multi-name)
- [ ] Maintain known-constants map (configurable per recipe)
- [ ] For known constants: remove PARAMETER + type decl, add to USE..ONLY
- [ ] For unknown constants: convert to declaration + executable assignment
- [ ] Detect and order PARAMETER dependency chains
- [ ] Flag PARAMETERs used in constant-expression contexts
- [ ] Remove old `DOUBLE PRECISION` declarations for converted PARAMETERs

### Phase 6: DATA Statement Conversion

**Goal:** Decompose DATA statements with FP values into declaration + executable
assignment, with explicit SAVE.

**Files affected:**
- `pyengine/fortran_migrator.py` — new function `convert_data_stmts()`
- `pyengine/flang_parser.py` — DATA context facts (from Phase 1.5)

**Implementation notes:**
- Variables initialized via DATA have implicit `SAVE` attribute. When converting
  to executable assignments, the migrator MUST add explicit `SAVE` to preserve
  semantics. Without this, variables would be re-initialized on every call —
  a silent correctness bug.
- Mixed DATA statements (some vars FP, some integer): only extract the FP parts,
  leave integer parts in DATA.
- Multi-line DATA statements: must join continuation lines, parse the full
  statement, then re-emit.
- Array initialization: `DATA (A(I), I=1,3) / 1.0, 2.0, 3.0 /` may need array
  constructor or explicit loop.
- Insert assignments at the start of the executable section (after all
  declarations, before the first executable statement).

**Inventory of DATA statements in BLAS/LAPACK:**
- `drotm.f` / `srotm.f`: `DATA ZERO,TWO/0.D0,2.D0/` — simple scalar FP
- `drotmg.f` / `srotmg.f`: `DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/` and
  `DATA GAM,GAMSQ,RGAMSQ/...` — scalar FP with non-trivial values
- `dlaln2.f` / `slaln2.f`: logical and integer DATA only (no FP migration needed)
- `dlaruv.f` / `slaruv.f`: integer array DATA only (no FP migration needed)
- `dlasy2.f` / `slasy2.f`: mixed — needs inspection

**Tasks:**
- [ ] Parse DATA statements (join continuation lines)
- [ ] Identify FP variables/values vs. non-FP (leave non-FP in DATA)
- [ ] Generate `TYPE(float64x2) :: name` declarations with explicit `SAVE`
- [ ] Generate assignment statements with constructor calls
- [ ] Insert at executable section boundary
- [ ] Handle array DATA initialization
- [ ] Handle mixed DATA statements (extract only FP parts)

### Phase 7: Intrinsic Function Mapping

**Goal:** Map intrinsic function calls to `multifloats` equivalents.

**Files affected:**
- `pyengine/intrinsics.py` — new mapping for multifloats mode
- `pyengine/fortran_migrator.py` — `replace_intrinsic_calls_multifloats()`

**Implementation notes:**
- The existing `INTRINSIC_MAP` has two modes: `needs_kind=False` (rename) and
  `needs_kind=True` (add KIND argument). Multifloats needs a third mode:
  "replace with constructor call" for conversion intrinsics.

  | Current intrinsic | KIND-based output | Multifloats output |
  |---|---|---|
  | `DBLE(x)` | `REAL(x, KIND=16)` | `float64x2(x)` |
  | `DCMPLX(a, b)` | `CMPLX(a, b, KIND=16)` | `complex128x2(a, b)` |
  | `DCOS(x)` | `COS(x)` | `COS(x)` (module generic) |
  | `DABS(x)` | `ABS(x)` | `ABS(x)` (module generic) |
  | `REAL(x)` (in expr) | `REAL(x, KIND=16)` | `float64x2(x)` |
  | `DIMAG(z)` | `AIMAG(z)` | `AIMAG(z)` (module generic) |

- The `multifloats` module must provide same-name generics (ABS, SQRT, SIN, COS,
  etc.) that accept `float64x2` arguments. If so, `USE multifloats` is sufficient
  and the simple-rename path works unchanged.
- For conversion intrinsics, implement a `needs_constructor=True` mode that emits
  `float64x2(inner)` instead of `REAL(inner, KIND=k)`.
- Update INTRINSIC declaration removal (these are no longer intrinsic when provided
  by a module).

**Tasks:**
- [ ] Define multifloats intrinsic mapping (third mode: constructor)
- [ ] Verify multifloats module provides same-name math generics
- [ ] Handle conversion intrinsics: DBLE, DCMPLX, REAL (in expression context)
- [ ] Handle `replace_generic_conversions()` equivalent for multifloats
- [ ] Update INTRINSIC declaration handling

### Phase 8: Free-Form and `la_constants` Support

**Goal:** Handle LAPACK's free-form `.f90` files and `la_constants` module
interaction.

**Files affected:**
- `pyengine/fortran_migrator.py` — free-form migration path

**Implementation notes:**
- Free-form LAPACK files use `integer, parameter :: wp = kind(1.d0)` and then
  `real(wp)` for declarations. For multifloats, the `wp` pattern is invalid —
  `TYPE(float64x2)` is not parameterized by an integer kind.
- Strategy: remove the `wp` parameter definition, replace `real(wp)` with
  `TYPE(float64x2)`, and update literal suffixes accordingly.
- `la_constants` rewriting (`rewrite_la_constants_use()`) currently maps to
  `LA_CONSTANTS_EP`. For multifloats, it should map to a multifloats-compatible
  constants module or import constants directly from the multifloats module.

**Tasks:**
- [ ] Handle `wp` parameter removal in free-form files
- [ ] Replace `real(wp)` / `complex(wp)` with `TYPE(float64x2)` / `TYPE(complex128x2)`
- [ ] Update `rewrite_la_constants_use()` for multifloats mode
- [ ] Update literal replacement for free-form (no column constraints, but
      constructor syntax still applies)

### Phase 9: Convergence Pipeline Update

**Goal:** The convergence/divergence checking and verification pipelines work
correctly for multifloats output.

**Files affected:**
- `pyengine/pipeline.py` — normalization and canonicalization functions

**Implementation notes:**
- `_canonicalize_for_compare()` must handle `TYPE(float64x2)` declarations and
  `float64x2(...)` constructor calls.
- `_light_normalize()` must handle `USE multifloats` statements.
- `cmd_verify()` residual checks must look for leftover `DOUBLE PRECISION`,
  D-exponent literals, and unconverted PARAMETER/DATA statements.
- Convergence normalization: `float64x2('1.0')` from a D-source and
  `float64x2('1.0')` from an S-source should be considered equivalent.

**Tasks:**
- [ ] Update `_canonicalize_for_compare()` for multifloats syntax
- [ ] Update `_light_normalize()` for multifloats declarations
- [ ] Update `cmd_verify()` residual checks
- [ ] Verify S/D and C/Z convergence under multifloats mode

### Phase 10: Build System and Integration Testing

**Goal:** End-to-end migration of BLAS, then LAPACK, with compilation.

**Files affected:**
- `pyengine/__main__.py` — CMake generation for multifloats

**Tasks:**
- [ ] Update CMake generation to link multifloats module
- [ ] Migrate all BLAS routines, verify compilation
- [ ] Run convergence checking under multifloats mode
- [ ] Run verification checks (no residual DOUBLE PRECISION)
- [ ] Extend to LAPACK
- [ ] Flag EQUIVALENCE blocks (`dlaln2.f`, `slaln2.f`) for manual review
- [ ] Document recipe configuration for multifloats target

---

## 6. Testing Strategy

### Unit Tests (per transformation function)

- **Type declarations:** All input forms (DOUBLE PRECISION, REAL*8, REAL(KIND=8),
  COMPLEX*16, DOUBLE COMPLEX), function return types, lines with attributes.
- **Literals:** Bare floats (`1.0`), D-exponent (`1.0D+0`), exact values,
  non-exact values, literals inside expressions, literals inside string literals
  (must be untouched), continuation lines.
- **PARAMETER conversion:** Known constants, unknown constants, chained PARAMETERs,
  constant-expression context (should warn).
- **DATA conversion:** Simple scalar, mixed FP/integer, SAVE attribute addition.
- **USE insertion:** Single subroutine, multiple subroutines per file,
  already-has-USE, fixed-form column constraints.
- **Intrinsics:** Conversion intrinsics (DBLE, DCMPLX), math intrinsics (DCOS),
  generic conversions (REAL(x) in expression context).

### Integration Tests (per file)

- `dgemv.f` -> validate against hand-migrated reference from Phase 0.
- `drotm.f` -> DATA conversion + SAVE.
- `dgemm.f` -> PARAMETER conversion.
- A free-form `.f90` file -> wp parameter path.

### End-to-End Tests

- Full BLAS migration with `--target multifloats`: compilation against stub module.
- Convergence check: S/D pairs and C/Z pairs converge.
- Verify check: no residual precision-specific keywords.

### Regression Tests

- KIND=10 and KIND=16 modes completely unaffected by multifloats additions.
- Existing convergence pipeline continues to pass for KIND-based migration.

---

## 7. Open Questions

1. **Routine naming convention:** Do multifloats-migrated routines use Q/X prefixes
   (same as KIND=16) or a new prefix (e.g., M/N)? This affects the rename map and
   library naming.

2. **`gfortran` parser support:** Can `gfortran_parser.py` extract the semantic
   information needed for multifloats, or is this flang-only?

3. **`multifloats` module API surface:** What operators, generics, and constants
   does the module actually provide? Phase 0 should define the minimum required API.

4. **Mixed-precision strategy:** When `float64x2` and `DOUBLE PRECISION` coexist
   in an expression (e.g., calling a DOUBLE PRECISION helper), what is the explicit
   conversion call? `DBLE(x)` for downcast? `float64x2(x)` for upcast?

5. **I/O statements:** Does the `multifloats` module provide defined I/O
   (`WRITE(unit, fmt) x` where x is `TYPE(float64x2)`), or must the migrator
   insert explicit conversions for output?

6. **Scope of initial target:** Should the first implementation target only
   fixed-form BLAS (simpler, well-understood), deferring LAPACK/free-form to a
   follow-up?
