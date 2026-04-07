# Review: `multifloats` Support Plan

## 1. Executive Summary

The proposed plan to support `multifloats` via a "Semantic Oracle + Source Rewriting" approach is architecturally sound and prioritizes the preservation of the original source's structure (macros, comments, formatting). This is critical for high-performance libraries like LAPACK. However, the implementation sequence needs re-ordering to integrate semantic information earlier, and several structural transformations (DATA/PARAMETER movement) require more robust handling of Fortran-specific edge cases like `SAVE` semantics and fixed-form constraints.

## 2. Architectural Critique

### Semantic Oracle vs. AST-to-AST
- **Strength:** The choice of "Option D" (Semantic Oracle) is superior to pure AST regeneration for this project. It ensures the output remains maintainable and "git-diff-able," which is essential for upstreaming.
- **Dependency:** The success of this approach hinges entirely on the quality of facts extracted from `flang-new`. If the oracle fails (e.g., due to missing module dependencies), the regex engine falls back to brittle heuristics.

### Phase Re-ordering
- **Critique:** Phase 8 (Semantic Oracle Enrichment) is scheduled too late.
- **Recommendation:** Move Semantic Oracle integration to **Phase 1**. The semantic facts are necessary to correctly implement Phase 2 (Type Declarations) and Phase 4 (Literal Replacement) without relying on fragile regex-based type inference.

## 3. Implementation Phases Critique

### Phase 0: Prototype Sufficiency
- **Observation:** `dgemv.f` is an excellent starting point but lacks `DATA` statements.
- **Recommendation:** Expand Phase 0 to include `drotm.f` or `dlarfg.f` to test the decomposition of `DATA` statements into executable assignments early in the prototype stage.

### Missing Phase: Module Mocking
- **Observation:** To run `flang-new` with sema on the *target* library, the compiler needs to see the `multifloats` module (or a stub).
- **Recommendation:** Add a phase for generating "Module Mocks" that allow the oracle to resolve the new types and generics even before the library is fully migrated.

## 4. Technical Risks & Mitigations

### Fixed-Form Line Length (72 Columns)
- **Risk:** `TYPE(float64x2)` is manageable (15 chars), but constructor calls like `float64x2('1.0D+0')` are significantly longer than `1.0D+0`. In fixed-form code, this will frequently trigger column-72 overflows.
- **Mitigation:** The rewriter must have a robust, column-aware line-breaking engine that can insert continuation characters (`+` or `&`) at safe points (e.g., after commas or operators) while respecting the 6-space indentation.

### SAVE Semantics in DATA Conversion
- **Risk:** In Fortran, variables initialized in `DATA` statements have an implicit `SAVE` attribute. Moving these to executable assignments removes this implicit attribute.
- **Mitigation:** The migrator must explicitly add the `SAVE` attribute to any variable whose `DATA` initialization is moved to the executable section.

### Mixed-Precision Promotion
- **Risk:** Standard Fortran promotes `INTEGER` to `REAL` automatically. Overloaded operators in `multifloats` must handle `TYPE(float64x2) + INTEGER`. However, `TYPE(float64x2) + DOUBLE PRECISION` might not be covered or might lead to precision loss if the compiler promotes the `multifloat` down to `double`.
- **Recommendation:** The migrator should identify mixed-precision expressions using the semantic oracle and insert explicit constructor calls (e.g., `A + float64x2(B)`) to ensure the desired precision is maintained.

## 5. Standards & Portability

### COMMON and EQUIVALENCE Blocks
- **Risk:** Derived types have strict restrictions in `COMMON` and `EQUIVALENCE` blocks. Specifically, a sequence-derived type is required for `COMMON`, and even then, mixing it with intrinsic types (like `INTEGER`) can be non-portable or prohibited.
- **Recommendation:** The migrator should flag `COMMON` or `EQUIVALENCE` blocks containing migrated variables. A manual review or a `BIND(C)` / `SEQUENCE` attribute in the `multifloats` module definition may be necessary.

### PARAMETER to Runtime Variables
- **Risk:** Converting `PARAMETER` to a runtime variable changes it from a compile-time constant to a mutable variable. This could break code that uses these "constants" in contexts requiring constant expressions (e.g., array bounds, other PARAMETER initializers).
- **Recommendation:** Use a two-tiered approach:
    1. Prefer `USE multifloats, ONLY: ZERO => MF_ZERO` for common constants.
    2. If a `PARAMETER` must be converted to a variable, verify it is not used in a constant-expression context.

## 6. Recommendations for Plan Update

1. **Prioritize Phase 8:** Integrate the Semantic Oracle facts immediately after Phase 1.
2. **Explicit SAVE:** Update Phase 6 to include the addition of the `SAVE` attribute.
3. **Module Mocking:** Add a task to create a minimal `.mod` or `.f90` stub for the oracle.
4. **Column-Aware Rewriter:** Ensure the `rewriter.cpp` is enhanced to handle fixed-form continuation for long constructor strings.
