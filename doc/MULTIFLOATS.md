# Revised Plan: `multifloats` Support (v3)

This document supersedes `04-plan-claude.md`, incorporating feedback from
`05-review-gemini.md`.

---

## 1. Key Changes from v2

1. **Simplified line-breaking algorithm** — the split mask now protects only
   string literal interiors, not all parenthesized expressions. The previous
   algorithm would have prevented splitting long BLAS call argument lists,
   causing unrecoverable column-72 overflows.
2. **Pattern B USE statement** — structural replacement via semantic oracle
   line-range identification, rather than fragile line-by-line regex surgery.
3. **Factory function exclusivity** — `TargetMode` should only be instantiated
   via factory functions (`kind_target`, `multifloats_target`) to enforce
   field consistency.
4. **DD/ZZ prefix convention** — multifloats routines use DD (real) and ZZ (complex)
   prefixes, distinct from Q/X (KIND=16).
5. **All open questions resolved** — gfortran parser support, external module
   availability, mixed-precision strategy, I/O conversion, initial scope.

---

## 2. `TargetMode` Design

### Consumption Patterns in the Codebase

Analysis of all 31 call sites where `kind: int` is consumed reveals six distinct
usage patterns:

| Pattern | Example | Sites |
|---|---|---|
| Type declaration output | `f'REAL(KIND={kind})'` | ~8 |
| Literal suffix | `f'{lit}E{exp}_{kind}'` | ~4 |
| Intrinsic KIND argument | `f'{name}({inner}, KIND={kind})'` | ~4 |
| Prefix map lookup | `PREFIX_MAP[kind]` | 4 |
| Kind parameter definition | `_KIND_PARAM_RE.sub(rf'\g<1>{kind}', line)` | 1 |
| la_constants rename map | `_LA_FULL_MAP[kind]` | 1 |

### Recommended Dataclass

```python
@dataclass(frozen=True)
class TargetMode:
    # Identity
    name: str                          # 'kind10', 'kind16', 'multifloats'

    # Type declaration output (Pattern 1)
    real_type: str                     # 'REAL(KIND=16)' or 'TYPE(float64x2)'
    complex_type: str                  # 'COMPLEX(KIND=16)' or 'TYPE(complex128x2)'

    # Literal replacement strategy (Pattern 2)
    literal_mode: str                  # 'kind_suffix' or 'constructor'
    kind_suffix: int | None            # 16 for KIND, None for multifloats
    real_constructor: str | None       # None for KIND, 'float64x2' for MF
    complex_constructor: str | None    # None for KIND, 'complex128x2' for MF

    # Intrinsic call transformation (Pattern 3)
    intrinsic_mode: str                # 'add_kind' or 'wrap_constructor'

    # Prefix map for renaming (Pattern 4)
    prefix_map: dict[str, str]         # {'R': 'W', 'C': 'U'}

    # Free-form kind parameter (Pattern 5)
    kind_param_value: str | None       # '16' for KIND, None for MF

    # Module / constants support (Pattern 6)
    module_name: str | None            # None for KIND, 'multifloats' for MF
    known_constants: dict[str, str]    # {'ZERO': 'MF_ZERO', ...}
    la_constants_map: dict[str, str]   # la_constants rename map

    @property
    def is_kind_based(self) -> bool:
        return self.kind_suffix is not None
```

**Why `frozen=True`:** TargetMode is a value object — immutable once constructed,
safely hashable for caching (e.g., the `_RENAME_PATTERN_CACHE`).

**Why not a class hierarchy:** The consumer code (transform functions) must make
per-call-site decisions about output format. A class hierarchy would push format
decisions into virtual methods, distributing transformation logic across subclasses
instead of keeping it in the transform functions where it belongs. A frozen
dataclass with factory functions is transparent and extensible.

### Factory Functions

```python
# In pyengine/target_mode.py
def kind_target(kind: int) -> TargetMode:
    """Construct TargetMode for KIND=10 or KIND=16."""
    ...

def multifloats_target(
    module: str = 'multifloats',
    real_type: str = 'float64x2',
    complex_type: str = 'complex128x2',
    known_constants: dict[str, str] | None = None,
    prefix_style: str = 'wide',   # DD/ZZ prefixes for float64x2
) -> TargetMode:
    """Construct TargetMode for multifloats."""
    ...
```

**Instantiation policy:** `TargetMode` should only be constructed through these
factory functions, never by direct `TargetMode(...)` calls. This ensures field
consistency (e.g., `literal_mode='constructor'` always paired with non-None
`real_constructor`). Enforce by convention (private `__init__` is not idiomatic
for frozen dataclasses, but all call sites should use the factories).

### Prefix Convention

**Multifloats uses DD/ZZ prefixes** — DD for real (`float64x2`), ZZ for complex
(`complex128x2`). The two-letter prefix avoids single-letter collisions and
makes the double-double nature explicit. These prefixes are distinct from all
existing conventions:

| Real | Complex | Precision |
|------|---------|-----------|
| S | C | Single (32-bit) |
| D | Z | Double (64-bit) |
| E | Y | Extended (80-bit, KIND=10) |
| Q | X | Quad (128-bit, KIND=16) |
| **DD** | **ZZ** | **Double-double (float64x2, multifloats)** |

Examples: `DDGEMM`, `DDAXPY`, `DDNRM2` (real); `ZZGEMM`, `ZZAXPY` (complex).

`PREFIX_MAP` remains keyed by int for backward compatibility. The `TargetMode`
carries its own `prefix_map` field, set by the factory function. The four
`PREFIX_MAP[kind]` lookup sites change to `target.prefix_map`.

### Migration Path (Three Sub-Phases)

**Phase 1a:** Introduce `TargetMode` dataclass and factory functions. Ship
the dataclass without changing any function signatures. (`1 PR`)

**Phase 1b:** Thread `TargetMode` through the top-level pipeline. CLI constructs
`kind_target(args.kind)` at the entry point. Pipeline functions accept
`TargetMode`. Inside pipeline functions, extract `target.kind_suffix` to pass to
leaf transform functions — an adapter layer. Add `--target multifloats` CLI
option. (`1 PR`)

**Phase 1c:** Update leaf transform functions one at a time to accept
`TargetMode` directly and branch on `target.is_kind_based`. Each function is a
separate, reviewable PR. At no point does KIND=10/16 behavior regress.

---

## 3. SAVE Semantics — Corrected Understanding

The previous plan (02) conflated two distinct Fortran behaviors:

### DATA + implicit SAVE
```fortran
      DATA X / 1.0D0 /
```
Initialization happens **exactly once** (before the first call). SAVE preserves
the value across calls. If X is later modified (`X = 2.0D0`), the modified value
persists to the next call.

### Executable Assignment + SAVE
```fortran
      TYPE(float64x2), SAVE :: X
      X = float64x2('1.0')   ! runs EVERY call
```
The assignment re-executes on **every call**. SAVE prevents deallocation between
calls but does NOT prevent re-execution. If X was modified, the assignment resets
it.

### Correct Migration Strategy

**Case 1: Variable is never written after initialization** (ZERO, ONE, TWO, GAM,
GAMSQ, RGAMSQ in BLAS). Re-initialization is harmless and idempotent. Either
approach works:
```fortran
      TYPE(float64x2) :: ZERO
      ZERO = float64x2(0.0D0)    ! re-executed each call, but idempotent
```
Adding SAVE is optional but semantically cleaner.

**Case 2: Variable IS written after initialization.** The DATA SAVE semantics
preserve the modified value across calls. An executable assignment would reset it.
This case requires a first-call guard or module-level initialization. However,
**no BLAS/LAPACK DATA-initialized FP variables are written after initialization**
— they are all used as read-only constants.

**Diagnostic:** The migrator should check whether any DATA-initialized FP variable
appears on the left-hand side of an assignment in the routine body. If so, emit a
warning for manual review.

---

## 4. Syntax-Aware Line Breaking

### Problem

`reformat_fixed_line()` at `fortran_migrator.py:556` does naive 66-character
chunking with backward scan for commas/spaces. A pathological case:

```fortran
      ABCDEFGHIJKLMNOPQRSTUVWXYZ1234=float64x2('1.0D+0')+float64x2('2
     +.0D+0')+ZZZZZ
```

This breaks the string literal `'2.0D+0'` across two lines — invalid Fortran.

### Solution

Add a `_build_split_mask()` helper that marks positions inside string literals
as unsafe for splitting:

```python
def _build_split_mask(body: str) -> list[bool]:
    """Return True where splitting is safe (outside string literals)."""
    mask = [True] * len(body)
    in_string = False
    quote_char = ''
    i = 0
    while i < len(body):
        ch = body[i]
        if in_string:
            mask[i] = False
            if ch == quote_char:
                if i + 1 < len(body) and body[i + 1] == quote_char:
                    mask[i + 1] = False  # doubled quote escape
                    i += 2
                    continue
                in_string = False
        elif ch in ("'", '"'):
            in_string = True
            quote_char = ch
            mask[i] = False
        i += 1
    return mask
```

**Design note:** The mask protects only string literal interiors, NOT all
parenthesized expressions. The v2 algorithm marked everything inside parentheses
as no-split, which would have prevented splitting long BLAS call argument lists
like:
```fortran
      CALL WGEMV('N', M, N, float64x2('1.0D+0'), A, LDA, X, 1, float64x2('0.0D+0'), Y, 1)
```
That line is 89 characters with commas between arguments that are perfectly valid
split points. The real risk is only splitting inside `'1.0D+0'` — which the
string-only mask prevents. Commas and spaces between function arguments remain
valid split points at any nesting depth. The existing backward scan for
commas/spaces naturally avoids splitting inside atomic tokens like `float64x2`
(no comma/space within the identifier).

Modify the split loop to check `safe[i]` before accepting a break point. This is
~20 lines of new code. Both call sites (`migrate_fixed_form` line 636 and
`_migrate_fixed_form_flang` line 968) benefit automatically.

### Free-Form Line Length

The free-form path has no line-length enforcement. Free-form Fortran has a
132-column limit. For multifloats, lines with multiple literal replacements could
conceivably exceed this. A defensive `reformat_free_line()` using `&` continuation
should be added, lower priority than the fixed-form fix.

---

## 5. Free-Form Migration: Two Distinct Patterns

### Pattern A: Self-Contained `wp` (BLAS `.f90` files)

Example: `drotg.f90`, `dnrm2.f90`

```fortran
integer, parameter :: wp = kind(1.d0)
real(wp) :: c, f, g, r, s
real(wp), parameter :: zero = 0.0_wp, one = 1.0_wp
real(wp), parameter :: safmin = real(radix(0._wp),wp)**max(...)
0.5_wp * x
```

**Multifloats transformation:**
1. Delete `integer, parameter :: wp = kind(1.d0)` entirely.
2. Insert `USE multifloats`.
3. Replace `real(wp)` -> `TYPE(float64x2)`, `complex(wp)` -> `TYPE(complex128x2)`.
4. Replace `_wp` literal suffixes: `0.5_wp` -> `float64x2(0.5D0)` or
   `float64x2('0.5D0')`.
5. Known PARAMETER constants (`zero`, `one`): import from module.
6. **Scaling PARAMETERs** (`safmin`, `safmax`, `tsml`, `tbig`): these use inquiry
   intrinsics (`radix(0._wp)`, `minexponent(0._wp)`, etc.) that are invalid on
   derived types. The entire scaling PARAMETER block must be replaced with module
   imports from a multifloats-compatible constants module.

### Pattern B: `la_constants`-Imported `wp` (LAPACK `.f90` files)

Example: `dlartg.f90`, `dlassq.f90`

```fortran
use LA_CONSTANTS, &
   only: wp=>dp, zero=>dzero, half=>dhalf, one=>done, &
         safmin=>dsafmin, safmax=>dsafmax
real(wp) :: x
```

**Multifloats transformation:**
1. Rewrite `USE LA_CONSTANTS` -> `USE la_constants_mf` (or equivalent).
2. Remove `wp=>dp` from the ONLY list (wp is no longer meaningful).
3. Rename remaining constants to multifloats equivalents.
4. Insert `USE multifloats` for type definitions.
5. Replace `real(wp)` -> `TYPE(float64x2)`.
6. Replace `_wp` literal suffixes as in Pattern A.

**Implementation strategy (revised):** The v2 plan noted that surgical removal
of `wp=>dp` from multi-line USE statements is difficult with line-by-line regex
when `wp=>dp` may be on a different continuation line than the module name.

**New approach:** Use semantic oracle facts to identify the exact line range of
the `USE LA_CONSTANTS` statement (including `&` continuations). Perform a
structural replacement of the entire USE block as a unit — parse the ONLY list,
remove the `wp=>dp` entry, rename the module and remaining constants, and re-emit
the block. This avoids fragile line-by-line surgery entirely.

### New Literal Form: `_wp` Suffixes

The existing `replace_literals()` handles D/E exponent forms and bare floats but
NOT `_wp` suffixed literals (`0.0_wp`, `1.0_wp`, `0.5_wp`). A new regex pattern
is needed:

```
(\d+\.\d*|\d*\.\d+)_wp  ->  float64x2(\1D0)   or   float64x2('\1D0')
```

This is specific to the free-form migration path and only applies in multifloats
mode (for KIND-based migration, changing `wp`'s value suffices).

### Inquiry Intrinsics

Patterns like `radix(0._wp)`, `minexponent(0._wp)`, `maxexponent(0._wp)`,
`digits(0._wp)`, `epsilon(0._wp)`, `huge(0.0_wp)` appear in PARAMETER
definitions for scaling constants. These intrinsics accept a value argument whose
type determines the result. They are not defined for derived types.

**Strategy:** The multifloats constants module must export pre-computed scaling
constants (`MF_SAFMIN`, `MF_SAFMAX`, etc.). The migrator replaces the entire
PARAMETER block with module imports.

### `LA_XISNAN` Module

`dlassq.f90` and others import `LA_XISNAN`. This needs a multifloats-compatible
equivalent (either the module handles `float64x2` via overloading, or the migrator
maps to a multifloats-specific variant).

---

## 6. PARAMETER Conversion — Refined

### Mixed-Type PARAMETER Lines

`dlaruv.f` line 115:
```fortran
      PARAMETER          ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )
```

Three parameters on one line: two INTEGER, one DOUBLE PRECISION. The migrator
must **partially decompose** this statement:
1. Identify FP entries via type declarations or semantic oracle.
2. Extract only FP entries (`R`).
3. Rewrite the PARAMETER line to contain only integer entries:
   `PARAMETER ( LV = 128, IPW2 = 4096 )`.
4. Add `R` as a runtime-initialized variable.

### PARAMETER Chaining

Dependency graph for `dlaruv.f`:
```
ONE = 1.0D0          (literal — import from module)
LV  = 128            (integer — keep as PARAMETER)
IPW2 = 4096          (integer — keep as PARAMETER)
R   = ONE / IPW2     (depends on ONE — runtime init)
```

`ONE` comes from module import (available immediately). `R = ONE / IPW2` goes
in the executable section before first use. Since `ONE / IPW2` is a cheap,
idempotent expression, re-computation on every call is acceptable.

### Complete PARAMETER Inventory (BLAS/LAPACK)

| File | Parameters | Type | Strategy |
|---|---|---|---|
| `dgemm.f` et al. | ONE, ZERO | Known constants | Module import |
| `dlaln2.f` | ZERO, ONE, TWO | Known constants | Module import |
| `dlasy2.f` | ZERO, ONE, TWO, HALF, EIGHT | Known + non-standard | HALF/EIGHT: runtime init or add to module |
| `dlaruv.f` | ONE, R (+ integers LV, IPW2) | Mixed-type, chained | Module import (ONE), runtime init (R), keep integers |
| Free-form files | zero, one, half, safmin, safmax, tsml, tbig, ssml, sbig | Known + scaling | Module import for all |

---

## 7. DATA Conversion — Refined

### Inventory with Semantic Analysis

| File | DATA Statement | Variables | FP? | Written After Init? | Migration |
|---|---|---|---|---|---|
| `drotm.f` | `DATA ZERO,TWO/0.D0,2.D0/` | ZERO, TWO | Yes | No | Executable assignment |
| `drotmg.f` | `DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/` | ZERO, ONE, TWO | Yes | No | Executable assignment |
| `drotmg.f` | `DATA GAM,GAMSQ,RGAMSQ/4096.D0,...,5.9604645D-8/` | GAM, GAMSQ, RGAMSQ | Yes | No | Executable assignment; RGAMSQ needs string constructor |
| `dlaln2.f` | `DATA ZSWAP/.../, RSWAP/.../, IPIVOT/.../` | ZSWAP, RSWAP, IPIVOT | No | — | Leave untouched |
| `dlaruv.f` | `DATA (MM(I,J),...) / ... /` | MM | No (integer) | — | Leave untouched |
| `dlasy2.f` | `DATA LOCU12/.../, XSWPIV/.../, BSWPIV/.../` | LOCU12, etc. | No | — | Leave untouched |

**Key finding:** In BLAS/LAPACK, DATA statements needing migration are limited to
`drotm.f`/`srotm.f` and `drotmg.f`/`srotmg.f`. All other DATA statements contain
only integer/logical values.

### RGAMSQ: Non-Exact Value

`drotmg.f` line 114:
```fortran
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
```

- `GAM = 4096.D0` — exact in double (2^12). Use `float64x2(4096.0D0)`.
- `GAMSQ = 16777216.D0` — exact in double (2^24). Use `float64x2(16777216.0D0)`.
- `RGAMSQ = 5.9604645D-8` — NOT exact in double. Must use string constructor:
  `float64x2('5.9604645D-8')`.

### `GAM**2` Pattern

`drotmg.f` uses `DD1 = DD1*GAM**2`. If GAM is `TYPE(float64x2)`, the module must
provide `operator(**)` for `(float64x2, INTEGER)`. This is standard for numeric
derived types but must be part of the module API requirements.

---

## 8. EQUIVALENCE — Mandatory Manual Intervention

### Affected Files

Only `dlaln2.f` and `slaln2.f` (the D/S pair).

```fortran
      DOUBLE PRECISION   CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
      EQUIVALENCE        ( CI( 1, 1 ), CIV( 1 ) ),
     $                   ( CR( 1, 1 ), CRV( 1 ) )
```

All four arrays are DOUBLE PRECISION and would be migrated to `TYPE(float64x2)`.
Fortran prohibits EQUIVALENCE with non-SEQUENCE derived types.

### Recommended Manual Rewrite

Eliminate EQUIVALENCE by replacing flat-index access `CRV(k)` with 2D indexing
using column-major mapping:
- `CRV(1)` -> `CR(1,1)`, `CRV(2)` -> `CR(2,1)`, `CRV(3)` -> `CR(1,2)`,
  `CRV(4)` -> `CR(2,2)`

Or add a small helper: `CR_FLAT(k) = CR(MOD(k-1,2)+1, (k-1)/2+1)`.

### Migrator Behavior

The migrator should:
1. Detect EQUIVALENCE statements containing FP variables (via semantic oracle
   or type declarations).
2. Emit a clear diagnostic: "EQUIVALENCE in dlaln2.f contains floating-point
   variables that cannot be migrated to derived types. Manual rewrite required."
3. Skip migration of the affected file (or migrate everything except EQUIVALENCE
   and flag the result as incomplete).

---

## 9. Convergence Pipeline Updates

### `_canonicalize_for_compare()` — Needs Significant Update

| Current Step | Change for Multifloats |
|---|---|
| D/E exponent normalization | Still needed |
| Strip `_N` kind suffixes | Also strip `float64x2(...)` constructor wrappers |
| Strip `REAL(expr, KIND=N)` casts | Also strip `float64x2(expr)` constructor-as-cast |
| Strip `(KIND=N)` on type specs | Also canonicalize `TYPE(float64x2)` to bare form |
| S/D/C/Z prefix collapse to `@` | No change |

### `_light_normalize()` — Needs Moderate Update

- Type declaration regex `r'REAL\(KIND=\d+\)'` must also match `TYPE(float64x2)`.
- `USE multifloats` statements should be stripped for convergence comparison.

### `cmd_verify()` — Needs Update

For multifloats mode, check for:
- Residual `DOUBLE PRECISION` not converted to `TYPE(float64x2)`.
- Residual D-exponent literals NOT wrapped in a constructor (bare `1.0D0` is a
  bug; `float64x2('1.0D0')` is correct).
- Residual `REAL*8`, `COMPLEX*16` not converted.
- Unconverted PARAMETER/DATA statements containing FP values.
- Column overflow (unchanged).

---

## 10. Multifloats Module — Minimum Required API

Derived from analysis of actual BLAS/LAPACK code patterns. The module is
implemented externally; Phase 0 validates this API against the actual module.

### Constructors
- `float64x2(double_precision_value)` — from exact double
- `float64x2('string_literal')` — from string (lossless)
- `float64x2(integer_value)` — from integer (needed for `DBLE(IT1)` patterns)
- `complex128x2(a, b)` — from two `float64x2` or two doubles

### Arithmetic Operators
- `+`, `-`, `*`, `/` — binary, both `(float64x2, float64x2)` and mixed
  `(float64x2, DOUBLE PRECISION)`, `(float64x2, INTEGER)`
- `**` with integer exponent — needed for `GAM**2` pattern
- Unary `-`

### Comparison Operators
- `==`, `/=`, `<`, `>`, `<=`, `>=` — both `(float64x2, float64x2)` and mixed
  `(float64x2, DOUBLE PRECISION)` (needed for `X(I).EQ.1.0D0` in `dlaruv.f`)

### Math Generics (same-name overloads)
- `ABS`, `SQRT`, `SIN`, `COS`, `TAN`, `EXP`, `LOG`, `LOG10`
- `AIMAG`, `CONJG`, `REAL` (extract real part of complex)
- `MIN`, `MAX`, `SIGN`
- `MOD`

### Conversion
- `DBLE(float64x2)` -> `DOUBLE PRECISION` (downcast)
- `mf_to_double(float64x2)` -> `DOUBLE PRECISION` (explicit I/O conversion;
  named distinctly so call sites can be found and removed once defined I/O is
  supported by the module)
- Defined assignment from `DOUBLE PRECISION` (implicit upcast)

### Named Constants
- `MF_ZERO`, `MF_ONE`, `MF_TWO`, `MF_HALF` (minimum for BLAS)
- `MF_EIGHT` (for `dlasy2.f`)
- Scaling constants: `MF_SAFMIN`, `MF_SAFMAX`, `MF_TSML`, `MF_TBIG`,
  `MF_SSML`, `MF_SBIG` (for free-form LAPACK)

---

## 11. Implementation Plan

### Phase 0: Manual Prototype

**Goal:** Validate the transformation catalog end-to-end.

**Targets:**
- `dgemv.f` — type decls, literals, intrinsics, USE insertion
- `drotm.f` — DATA conversion (ZERO, TWO)
- `drotmg.f` — DATA with non-exact values (RGAMSQ), `INTRINSIC DABS`, `GAM**2`
- `dgemm.f` — PARAMETER conversion (ONE, ZERO)
- `drotg.f90` — free-form Pattern A (self-contained wp, scaling PARAMETERs)

**Deliverables:**
- [ ] Hand-migrated files for each target
- [ ] Transformation catalog (before/after for every change)
- [ ] Validate API (Section 10) against the external `multifloats` module
- [ ] Verification: migrated files compile against the module
- [ ] Gap analysis: any transformations not covered in this plan

### Phase 1: Target Mode Infrastructure

**Goal:** The migrator accepts `--target multifloats` without breaking KIND=10/16.

**Sub-phase 1a** (1 PR):
- [ ] Create `pyengine/target_mode.py` with `TargetMode` dataclass
- [ ] Implement `kind_target(kind)` and `multifloats_target(...)` factories
- [ ] Unit tests for both factories

**Sub-phase 1b** (1 PR):
- [ ] Update CLI: add `--target` option alongside `--kind`
- [ ] Construct `TargetMode` at CLI entry point via factory functions
- [ ] Thread `TargetMode` through pipeline functions
- [ ] Adapter layer: leaf functions still accept `kind: int` via `target.kind_suffix`
- [ ] `PREFIX_MAP` lookups -> `target.prefix_map`
- [ ] Regression tests: KIND=10/16 produce identical output

**Sub-phase 1c** (per-function PRs):
- [ ] Update `replace_type_decls()` to accept `TargetMode`, branch on `is_kind_based`
- [ ] Update `replace_literals()` — same pattern
- [ ] Update `replace_intrinsic_calls()` — same pattern
- [ ] Update `replace_generic_conversions()` — same pattern
- [ ] Update `_replace_kind_parameter()` — same pattern
- [ ] Update `rewrite_la_constants_use()` — same pattern
- [ ] Update `_build_sub_vars()` in `c_migrator.py` — same pattern

### Phase 1.5: Semantic Oracle Enrichment

**Goal:** Extract PARAMETER, DATA, and procedure boundary facts.

**Files affected:** `pyengine/flang_parser.py`, `pyengine/gfortran_parser.py`

**New `ParseTreeFacts` fields:**
- `parameter_stmts: list[ParameterInfo]` — name, value, type (FP/integer),
  line number
- `data_stmts: list[DataInfo]` — variable names, values, types, line numbers
- `procedure_boundaries: list[ProcInfo]` — name, header_line, first_exec_line,
  end_line
- `variable_types: dict[str, str]` — name -> type (for mixed-type disambiguation)
- `use_stmt_ranges: list[UseStmtInfo]` — module name, line range, ONLY list
  entries (for Pattern B structural replacement)

**Tasks:**
- [ ] Add `use_sema: bool` parameter to `run_flang_parse_tree()`
- [ ] When `use_sema=True`, invoke `-fdebug-dump-parse-tree` (with sema)
- [ ] Provide multifloats `.mod` via `-I` flag
- [ ] Extract PARAMETER/DATA/procedure/USE-statement facts
- [ ] Handle sema failures: fall back to no-sema + regex heuristics
- [ ] Implement equivalent extraction in `gfortran_parser.py`

### Phase 2: Type Declaration Rewriting

**Goal:** `DOUBLE PRECISION` -> `TYPE(float64x2)`, etc.

**Implementation:** New branch in `replace_type_decls()` activated when
`target.is_kind_based is False`. Emits `target.real_type` / `target.complex_type`.

**Tasks:**
- [ ] New regex patterns for `TYPE(...)` output syntax
- [ ] All input forms: `DOUBLE PRECISION`, `REAL*8`, `REAL(KIND=8)`, `REAL(8)`,
      `DOUBLE COMPLEX`, `COMPLEX*16`, `COMPLEX(KIND=8)`
- [ ] Function return types: `DOUBLE PRECISION FUNCTION FOO` ->
      `TYPE(float64x2) FUNCTION FOO`
- [ ] Preserve attributes (INTENT, DIMENSION, ALLOCATABLE, etc.)
- [ ] `replace_standalone_real_complex()` multifloats branch

### Phase 3: USE Statement Insertion

**Goal:** Every procedure using `float64x2` gets `USE multifloats`.

**Implementation:** New function `insert_use_multifloats(source, needed_constants)`
as a whole-source pass (not line-by-line).

**Tasks:**
- [ ] Detect procedure boundaries (SUBROUTINE/FUNCTION/PROGRAM)
- [ ] Insert `USE multifloats` after header, before IMPLICIT NONE
- [ ] Build ONLY clause dynamically (type names + needed constants)
- [ ] Avoid duplicate insertion
- [ ] Fixed-form column constraints (`      USE multifloats` = 22 chars, fits)

### Phase 4: Literal Replacement

**Goal:** `1.0D+0` -> `float64x2('1.0D+0')` or `float64x2(1.0D0)`.

**Implementation:** New branch in `replace_literals()` for constructor mode.

**Fixed-form literals:**
- D-exponent: `1.0D+0` -> `float64x2('1.0D+0')` (string) or
  `float64x2(1.0D0)` (exact)
- Bare floats: `1.0` -> `float64x2(1.0D0)` (exact)
- Must NOT replace inside PARAMETER/DATA statements (Phase 5/6 handle those)
- Must NOT replace integer/character literals

**Free-form literals** (additional pattern):
- `_wp` suffix: `0.5_wp` -> `float64x2(0.5D0)` or `float64x2('0.5D0')`
- New regex: `(\d+\.\d*|\d*\.\d+)_wp`

**Line breaking enhancement:**
- [ ] Implement `_build_split_mask()` — string-literal-only protection (Section 4)
- [ ] Modify `reformat_fixed_line()` to use the mask
- [ ] Verify no self-interference (constructor string literals not re-processed)

### Phase 5: PARAMETER Statement Conversion

**Goal:** Remove FP PARAMETER declarations; replace with module imports or
runtime initialization.

**Implementation:** New whole-source pass `convert_parameter_stmts()`.

**Tasks:**
- [ ] Parse PARAMETER statements (multi-line, multi-name)
- [ ] Mixed-type PARAMETER lines: extract only FP entries, preserve integers
- [ ] Known constants: remove PARAMETER + type decl, add to USE..ONLY
- [ ] Unknown constants: convert to `TYPE(float64x2) :: name` + exec assignment
- [ ] Chained PARAMETERs: topological sort, emit assignments in dependency order
- [ ] Scaling PARAMETERs (free-form `safmin`, `safmax`, etc.): replace entire
      block with module imports
- [ ] Diagnostic: flag PARAMETERs used in constant-expression contexts
- [ ] Handle `real(wp), parameter :: zero = 0.0_wp` in free-form (combined
      type-decl + PARAMETER + literal)

### Phase 6: DATA Statement Conversion

**Goal:** Decompose FP DATA statements into executable assignments.

**Implementation:** New whole-source pass `convert_data_stmts()`.

**Tasks:**
- [ ] Join continuation lines for multi-line DATA statements
- [ ] Identify FP variables via type declarations or semantic oracle
- [ ] Mixed DATA: extract only FP parts, leave integer/logical parts
- [ ] Generate declarations: `TYPE(float64x2) :: name`
- [ ] Generate executable assignments with constructor calls
- [ ] String constructor for non-exact values (RGAMSQ)
- [ ] Insert assignments at executable section boundary
- [ ] Diagnostic: warn if DATA-initialized variable is later written in routine

### Phase 7: Intrinsic Function Mapping

**Goal:** Map intrinsic calls to multifloats equivalents.

**Implementation:** Third intrinsic mode `needs_constructor=True`.

| Intrinsic | Mode | Output |
|---|---|---|
| `DBLE(x)` | constructor | `float64x2(x)` |
| `DCMPLX(a,b)` | constructor | `complex128x2(a, b)` |
| `DCOS(x)` | rename | `COS(x)` (module generic) |
| `DABS(x)` | rename | `ABS(x)` (module generic) |
| `REAL(x)` in expr | constructor | `float64x2(x)` |
| `REAL(z)` (complex->real part) | rename | `REAL(z)` (module generic) |
| `DIMAG(z)` | rename | `AIMAG(z)` (module generic) |

**Caution:** `REAL(t)` extracting the real part of a complex number (e.g., in
`zlartg.f90` statement function `ABSSQ(t) = real(t)**2 + aimag(t)**2`) is NOT a
type conversion. The migrator must distinguish real-part extraction from
type-conversion by context (argument type). The semantic oracle helps here.

**Tasks:**
- [ ] Add `needs_constructor` to `INTRINSIC_MAP` entries
- [ ] Implement constructor replacement in `replace_intrinsic_calls()`
- [ ] Handle `replace_generic_conversions()` for multifloats
- [ ] Remove `INTRINSIC` declarations for module-provided generics
- [ ] Distinguish `REAL(complex)` (real-part) from `REAL(integer)` (conversion)

### Phase 8: Free-Form and `la_constants` Support

**Goal:** Handle both free-form patterns (A and B) and la_constants interaction.

**Tasks:**
- [ ] Pattern A: Delete `wp` parameter line, insert USE, replace `real(wp)`,
      replace `_wp` suffixes, replace scaling PARAMETERs with module imports
- [ ] Pattern B: Use semantic oracle `use_stmt_ranges` to identify USE block
      line range; perform structural replacement of entire ONLY list (remove
      `wp=>dp`, rename module, rename constants) as a unit
- [ ] Handle `LA_XISNAN` -> multifloats equivalent
- [ ] Statement functions with `real(wp)` return type
- [ ] Function return type declarations (`real(wp) :: DNRM2`)
- [ ] Optional: add `reformat_free_line()` for 132-column protection

### Phase 9: Convergence Pipeline Update

**Tasks:**
- [ ] `_canonicalize_for_compare()`: strip `float64x2(...)` wrappers, canonicalize
      `TYPE(float64x2)` declarations
- [ ] `_light_normalize()`: match `TYPE(float64x2|complex128x2)` in type regexes,
      strip `USE multifloats`
- [ ] `cmd_verify()`: check for residual DOUBLE PRECISION, bare D-literals,
      unconverted PARAMETER/DATA, REAL*8/COMPLEX*16

### Phase 10: Build System and Integration Testing

**Tasks:**
- [ ] Update CMake generation: link multifloats module, add `-I` for `.mod` files
- [ ] Full BLAS migration: compile against external multifloats module
- [ ] Convergence check: S/D and C/Z pairs converge under multifloats
- [ ] Verify check: no residual precision-specific keywords
- [ ] Flag `dlaln2.f`/`slaln2.f` EQUIVALENCE for manual review
- [ ] Extend to LAPACK
- [ ] Document recipe configuration for multifloats target

---

## 12. Testing Strategy

### Unit Tests

Per transformation function, covering:
- All input declaration forms -> `TYPE(float64x2)` output
- D-exponent, bare-float, and `_wp`-suffix literal forms -> constructor output
- Exact-in-double classification (powers of 2, 0, 1, -1, 0.5)
- PARAMETER: known constants, unknown constants, chained, mixed-type lines
- DATA: simple scalar, non-exact values, SAVE semantics
- USE insertion: single/multiple procedures, duplicates, fixed-form constraints
- Intrinsics: all three modes (rename, add_kind, wrap_constructor)
- Line breaking: constructor calls at boundary, pathological cases, long BLAS
  calls with constructors in argument lists
- Self-interference: `float64x2('1.0D+0')` not double-processed

### Integration Tests

- `dgemv.f` -> validate against Phase 0 hand-migrated reference
- `drotm.f` -> DATA conversion
- `drotmg.f` -> non-exact DATA values, `GAM**2`, `INTRINSIC DABS`
- `dgemm.f` -> PARAMETER conversion
- `drotg.f90` -> free-form Pattern A, scaling constants
- `dlartg.f90` -> free-form Pattern B, la_constants, statement function

### End-to-End Tests

- Full BLAS migration with `--target multifloats`: compilation succeeds
- S/D and C/Z convergence passes
- Verify: no residual precision-specific keywords

### Regression Tests

- KIND=10 and KIND=16 produce identical output before and after all changes
- Existing convergence pipeline passes unchanged

---

## 13. Resolved Open Questions

1. **Routine naming convention:** Use DD/ZZ prefixes — DD for real (`float64x2`),
   ZZ for complex (`complex128x2`). The two-letter prefix avoids single-letter
   collisions. Distinct from Q/X (KIND=16) to avoid ambiguity between quad
   and double-double.

2. **EQUIVALENCE handling:** Mandatory manual intervention for `dlaln2.f`/
   `slaln2.f`. Concrete rewrite strategy: replace flat-index CRV(k) with 2D
   CR(i,j) using column-major mapping.

3. **SAVE semantics:** For BLAS/LAPACK, all DATA-initialized FP variables are
   read-only constants. Executable re-assignment is safe. Add diagnostic for
   variables written after initialization.

4. **`gfortran` parser support:** Yes — `gfortran_parser.py` should be able to
   extract PARAMETER/DATA context and variable type information. Both parser
   backends (flang and gfortran) should support multifloats mode.

5. **`multifloats` module availability:** The module is already implemented
   externally. Phase 0 validates the minimum API (Section 10) against the actual
   module rather than defining a new specification.

6. **Mixed-precision expressions:** Compare the single-precision and
   double-precision versions of each routine. Differences between the two
   indicate where type migration should occur. This is already the core strategy
   of the S/D and C/Z convergence pipeline — the same approach identifies which
   expressions involve precision-dependent types and thus need explicit conversion
   in multifloats mode.

7. **I/O statements:** The migrator inserts explicit conversion calls at I/O
   boundaries. These use a dedicated wrapper function `mf_to_double(x)` that
   converts `float64x2` to `DOUBLE PRECISION` for output:
   ```fortran
   ! Before
   WRITE(*,*) X

   ! After (multifloats)
   WRITE(*,*) mf_to_double(X)
   ```
   The `mf_to_double` name is chosen so that these call sites are easily
   searchable and can be mechanically removed once the multifloats module adds
   defined I/O support. The external library will provide `mf_to_double` as part
   of the module API.

8. **Scope of initial target:** Fixed-form BLAS first (Phases 0-7), then
   free-form LAPACK (Phase 8), then convergence/build (Phases 9-10).
