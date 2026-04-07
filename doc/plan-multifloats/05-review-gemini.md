# Review: Revised `multifloats` Support Plan v2 (`04-plan-claude.md`)

## 1. Executive Summary

The v2 plan for `multifloats` support is a highly mature and comprehensive document that addresses the deepest technical challenges of the migration. By identifying specific consumption patterns (31 call sites) and refining the `TargetMode` abstraction, the plan provides a clear, zero-regression path for modernizing the migrator's core infrastructure. The corrected understanding of `SAVE` semantics and the detailed strategies for free-form migration (Patterns A and B) significantly reduce the risk of silent correctness bugs.

## 2. Architectural Design Critique

### `TargetMode` and Migration Path
- **Strength:** The use of a **frozen dataclass** is the optimal architectural choice for immutability and reliability within the transformation caches.
- **Safety:** The **three-sub-phase migration path (1a, 1b, 1c)**, particularly the use of an adapter layer in Phase 1b, is an excellent strategy for ensuring that `KIND=10/16` functionality remains entirely unaffected during the transition.
- **Refinement:** The factory functions (`kind_target` and `multifloats_target`) should be the only way to instantiate these objects to maintain consistency across the codebase.

### Semantic Oracle (Phase 1.5)
- **Strength:** The plan correctly identifies that `flang-new` with full semantics is required to resolve ambiguities that regex cannot handle, such as mixed-type `PARAMETER` statements and `DATA` statements with non-FP variables.
- **Verification:** Ensure that the `-I` flag correctly points to the `multifloats.mod` (or stub) during the oracle's execution to resolve the new types.

## 3. Technical Risks & Mitigations

### Syntax-Aware Line Breaking (`_build_split_mask`)
- **Risk:** The proposed mask forbids **all** splitting inside parentheses (`depth > 0`). While this protects `float64x2(...)` constructors, it may prevent necessary splits in extremely long expressions or function calls, leading to column-72 overflows.
- **Refinement:** The mask should be refined to allow splitting at **commas or operators** even when `depth > 0`, provided the split point is not inside a nested string literal. The primary goal is to protect string literals and atomicity of tokens, not to forbid all intra-parenthetical splits.

### `SAVE` and `DATA` Semantics
- **Strength:** The distinction between read-only constants (where idempotent re-initialization is safe) and modified variables is a critical insight.
- **Verification:** The proposed diagnostic for manual review of variables appearing on the LHS of assignments is a necessary safeguard for numerical library migration.

### `EQUIVALENCE` and Mixed-Type `PARAMETER`
- **Strength:** Upgrading `EQUIVALENCE` to "Mandatory Manual Intervention" is the correct engineering decision. The partial decomposition strategy for `PARAMETER` lines is technically sound and necessary for the complex statements found in `dlaruv.f`.

## 4. Gaps and Refinements

### Free-Form Pattern B and `la_constants`
- **Observation:** Surgical removal of `wp=>dp` in Pattern B is complex due to multi-line `USE` statements.
- **Recommendation:** The migrator should use the semantic oracle facts to identify the exact line range of the `USE` statement and perform a structural replacement of the `ONLY` list rather than a simple line-by-line regex.

### `_wp` Literal Suffixes
- **Strength:** Adding a specific regex for `_wp` suffixes is essential for the free-form migration path (Pattern A), as these are pervasive in modern LAPACK routines.

## 5. Testing Strategy

- **Strength:** The multi-tiered testing strategy (Unit, Integration, E2E, Regression) is comprehensive.
- **Prototyping:** The expanded Phase 0 targets (`dgemv.f`, `drotm.f`, `drotmg.f`, `dgemm.f`, `drotg.f90`) correctly exercise the most difficult transformations, including `DATA` conversion and free-form scaling constants.

## 6. Conclusion

The v2 plan is exceptionally high quality and is **ready for implementation**. The most critical next steps are:
1.  **Phase 0 Execution:** Begin the manual migration of the expanded target set to validate the transformation catalog.
2.  **Phase 1a/b Implementation:** Introduce the `TargetMode` dataclass and thread it through the pipeline.
3.  **Line Breaking Refinement:** Update the split-mask logic to allow splits at commas/operators inside parentheses while protecting string literals.
