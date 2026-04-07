# Plan: `multifloats` Support (`float64x2` / `complex128x2`)

## Goal

Extend the migrator to target `float64x2` and `complex128x2` derived types
from an external `multifloats` module, in addition to the existing `KIND=10`
and `KIND=16` intrinsic-type targets.

## Key Differences from KIND-Based Migration

| Aspect | KIND=16 (current) | `multifloats` (new) |
|---|---|---|
| Type syntax | `REAL(KIND=16)` | `TYPE(float64x2)` |
| Literals | `1.0E+0_16` | `float64x2('1.0')` or `float64x2(1.0D0)` |
| Intrinsics | Generic resolution by compiler | Module-provided generics via `USE` |
| `PARAMETER` | Works (intrinsic constant expr) | Not a constant expression — must convert |
| `DATA` | Works | Derived types not allowed — must convert |
| Operators | Built-in | Overloaded via `USE multifloats` |
| Mixed precision | Compiler promotes automatically | Explicit conversion calls needed |

## Design Decision: Semantic Oracle + Source Rewriting (Option D)

**Do not modify parse trees and regenerate source.** Preprocessor directives
would be destroyed.

Instead, enrich the existing architecture:

```
Pass 1:  flang-new -fc1 -fdebug-dump-parse-tree  (WITH sema)
         ─► extract semantic facts:
            • variable types (which vars are DOUBLE PRECISION?)
            • expression types (for mixed-precision detection)
            • PARAMETER context (which literals are in PARAMETER stmts?)
            • DATA context (which literals are in DATA stmts?)
            • intrinsic vs. external call classification

Pass 2:  existing regex engine applies transformations
         to original source (macros intact)
         guided by semantic facts from Pass 1
```

## Literal Conversion Strategy

- **String constructor** `float64x2('1.0D+0')` as the safe default (no
  precision loss through intermediate representation).
- **Double constructor** `float64x2(1.0D0)` for values exactly representable
  in double precision (`0.0`, `1.0`, `2.0`, `-1.0`, etc.).
- The `multifloats` module provides both constructors.

## PARAMETER / DATA Conversion Strategy

### PARAMETER statements

`PARAMETER` with derived types is not valid standard Fortran (the initializer
must be a constant expression, and a derived-type constructor generally is not).

**Strategy:** The `multifloats` module exports common named constants:
`MF_ZERO`, `MF_ONE`, `MF_HALF`, `MF_TWO`, etc. (or `ZERO`, `ONE`, ... under
a renamed import).

Transformation:
```fortran
! Before
DOUBLE PRECISION ZERO, ONE
PARAMETER (ZERO=0.0D+0, ONE=1.0D+0)

! After
USE multifloats, ONLY: ZERO => MF_ZERO, ONE => MF_ONE
! (PARAMETER lines removed)
```

For non-standard constants not pre-defined in the module, convert to
runtime-initialized variables:
```fortran
! Before
DOUBLE PRECISION ALPHA
PARAMETER (ALPHA=0.9007199254740992D+15)

! After
TYPE(float64x2) :: ALPHA
ALPHA = float64x2('0.9007199254740992D+15')  ! in executable section
```

### DATA statements

`DATA` cannot initialize derived types.

Transformation:
```fortran
! Before
DOUBLE PRECISION X, Y
DATA X, Y / 1.0D0, 2.0D0 /

! After
TYPE(float64x2) :: X, Y
X = float64x2('1.0')   ! moved to executable section
Y = float64x2('2.0')   ! moved to executable section
! (DATA line removed)
```

**Implementation note:** This is a structural transformation — code moves from
the declaration section to the executable section. The migrator must track the
boundary between declarations and executable statements (the first executable
statement, or an implicit boundary after all declarations).

## Implementation Phases

### Phase 0: Prototype with `dgemv`

**Goal:** End-to-end manual migration of `dgemv.f` to validate the
transformation catalog.

Tasks:
- [ ] Hand-migrate `dgemv.f` to use `float64x2` / `multifloats`
- [ ] Catalog every transformation applied (type, line, before/after)
- [ ] Identify any transformation categories not covered below
- [ ] Create a `multifloats` stub module sufficient for `dgemv` to compile

### Phase 1: New Target Kind Infrastructure

**Goal:** The migrator accepts `--kind multifloats` (or similar) and activates
derived-type-aware transformation paths.

Files affected:
- `pyengine/__main__.py` — new CLI option / kind value
- `pyengine/pipeline.py` — dispatch to multifloats-aware transforms
- `pyengine/fortran_migrator.py` — new transformation mode

Tasks:
- [ ] Add `multifloats` as a recognized kind target
- [ ] Add configuration for multifloats module name, type names, and
      constructor names
- [ ] Gate existing KIND-based literal/type logic behind a kind-mode check

### Phase 2: Type Declaration Rewriting

**Goal:** `DOUBLE PRECISION` / `REAL(KIND=n)` → `TYPE(float64x2)`,
`COMPLEX*16` / `COMPLEX(KIND=n)` → `TYPE(complex128x2)`.

Files affected:
- `pyengine/fortran_migrator.py` — `replace_type_declarations()`

Tasks:
- [ ] New regex patterns for `TYPE(float64x2)` output (differs from
      `REAL(KIND=N)` syntactically)
- [ ] Handle all input forms: `DOUBLE PRECISION`, `REAL*8`, `REAL(KIND=8)`,
      `REAL(8)`, free-form `real(wp)`
- [ ] Handle function return type declarations:
      `DOUBLE PRECISION FUNCTION FOO` → `TYPE(float64x2) FUNCTION FOO`
- [ ] Preserve attributes: `INTENT`, `DIMENSION`, `ALLOCATABLE`, etc.

### Phase 3: USE Statement Insertion

**Goal:** Every procedure that uses `float64x2` gets `USE multifloats`.

Files affected:
- `pyengine/fortran_migrator.py` — new pass or hook in existing pipeline

Tasks:
- [ ] Detect procedure boundaries (SUBROUTINE/FUNCTION/PROGRAM)
- [ ] Insert `USE multifloats` after the procedure header, before other
      declarations
- [ ] Avoid duplicate insertion if already present
- [ ] Handle fixed-form column constraints for the inserted line
- [ ] Import needed constants (ZERO, ONE, etc.) via ONLY clause if
      PARAMETER lines were removed (coordinate with Phase 5)

### Phase 4: Literal Replacement

**Goal:** `1.0D+0` → `float64x2('1.0D+0')` (or double constructor for
exact values).

Files affected:
- `pyengine/fortran_migrator.py` — `replace_literals()`

Tasks:
- [ ] New literal replacement logic for derived-type mode
- [ ] Classify literals as exact-in-double vs. requires-string-constructor
- [ ] Handle literals inside expressions: `A = 1.0D0 * B` →
      `A = float64x2('1.0') * B` (operator overloading handles the rest)
- [ ] Do NOT replace literals inside PARAMETER statements (handled in Phase 5)
- [ ] Do NOT replace integer literals or character literals

### Phase 5: PARAMETER Statement Conversion

**Goal:** Remove `PARAMETER` declarations of floating-point constants and
replace with module imports or runtime initialization.

Files affected:
- `pyengine/fortran_migrator.py` — new pass
- `pyengine/flang_parser.py` — extract PARAMETER context facts

Tasks:
- [ ] Parser: identify which names are declared as PARAMETER with FP values
- [ ] Maintain a known-constants map (`ZERO` → `MF_ZERO`, `ONE` → `MF_ONE`,
      etc.) configurable per recipe
- [ ] For known constants: remove PARAMETER line, add to USE..ONLY import list
- [ ] For unknown constants: convert to `TYPE(float64x2) :: name` declaration
      + assignment in executable section
- [ ] Remove the old `DOUBLE PRECISION` declaration for converted PARAMETERs

### Phase 6: DATA Statement Conversion

**Goal:** Decompose DATA statements with floating-point values into
declaration + executable assignment.

Files affected:
- `pyengine/fortran_migrator.py` — new pass
- `pyengine/flang_parser.py` — extract DATA context facts

Tasks:
- [ ] Parser: identify DATA statements containing FP variables/values
- [ ] Parse DATA statement structure (variable lists / value lists)
- [ ] Handle mixed DATA statements (some vars FP, some integer) — only
      extract the FP parts, leave integer parts in DATA
- [ ] Generate assignment statements and insert at start of executable section
- [ ] Handle array initialization in DATA: `DATA (A(I), I=1,3) / 1.0, 2.0, 3.0 /`
      — may need array constructor or loop

### Phase 7: Intrinsic Function Mapping

**Goal:** Ensure intrinsic function calls resolve to `multifloats` generics.

Files affected:
- `pyengine/intrinsics.py` — new/modified map
- `pyengine/fortran_migrator.py` — `replace_intrinsic_calls()`

Tasks:
- [ ] Verify that `multifloats` provides same-name generics (ABS, SQRT, SIN,
      COS, etc.) — if so, `USE multifloats` is sufficient and current
      generic-ification logic works as-is
- [ ] For conversion intrinsics (`DBLE`, `REAL`, `DCMPLX`, `CMPLX`): replace
      with `float64x2(x)` / `complex128x2(x)` constructor calls instead of
      `REAL(x, KIND=16)`
- [ ] Handle `AIMAG`, `CONJG`, `DREAL`/`DIMAG` — must map to multifloats
      equivalents
- [ ] Update INTRINSIC declaration removal (these are no longer intrinsic)

### Phase 8: Semantic Oracle Enrichment

**Goal:** Switch from `-fdebug-dump-parse-tree-no-sema` to sema-enabled
parse tree when targeting multifloats, to extract type information.

Files affected:
- `pyengine/flang_parser.py` — parser invocation and fact extraction

Tasks:
- [ ] Enable sema in flang invocation for multifloats mode
- [ ] Extract variable type map (name → declared type)
- [ ] Extract expression type info where needed (for mixed-precision)
- [ ] Handle cases where sema fails (e.g., missing module dependencies) —
      fall back to no-sema + heuristics
- [ ] Ensure parse succeeds despite the source not yet being migrated
      (parser sees original DOUBLE PRECISION code)

### Phase 9: Integration, Testing, and Convergence

**Goal:** End-to-end migration of BLAS, then LAPACK.

Tasks:
- [ ] Migrate all BLAS routines, verify compilation with multifloats module
- [ ] Run convergence checking (S/C vs D/Z) under multifloats mode
- [ ] Run verification checks (no residual DOUBLE PRECISION)
- [ ] Extend to LAPACK — identify new transformation categories if any
- [ ] Document recipe configuration for multifloats target

## Risks and Open Questions

1. **Sema availability:** `flang-new` with sema requires that all referenced
   modules are available. For LAPACK this means building `LA_CONSTANTS` first.
   May need a bootstrapping order.

2. **Fixed-form line length:** `TYPE(float64x2)` is 16 chars vs. `REAL(KIND=16)`
   at 12 or `DOUBLE PRECISION` at 16. Constructor calls like
   `float64x2('1.0D+0')` are longer than `1.0E+0_16`. More continuation-line
   splitting will be needed.

3. **Mixed-precision expressions:** Where `DOUBLE PRECISION` and `float64x2`
   coexist in an expression (e.g., integer-to-float promotion, or calling a
   DOUBLE PRECISION helper), explicit conversion calls are needed. The semantic
   oracle helps identify these, but the transformation logic is non-trivial.

4. **Array constructors:** `(/ 1.0D0, 2.0D0 /)` must become
   `(/ float64x2('1.0'), float64x2('2.0') /)` — the array constructor must
   have a consistent type.

5. **I/O statements:** `WRITE(*,*) X` where X is `TYPE(float64x2)` needs
   defined I/O in the module, or explicit conversion to DOUBLE PRECISION for
   output. This is a module-design question, not a migrator question.

6. **Scope of `multifloats` module:** The module must provide a comprehensive
   set of operators, intrinsic-equivalent functions, and constants. The migrator
   assumes this exists; Phase 0 should validate the required API surface.
