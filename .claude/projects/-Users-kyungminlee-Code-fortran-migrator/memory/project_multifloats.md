---
name: multifloats module is external
description: The multifloats (float64x2/complex128x2) module is an externally implemented library, not part of this project
type: project
---

The `multifloats` module providing `float64x2` and `complex128x2` types is an externally implemented library. Phase 0 of the plan validates the migrator's minimum API requirements against the actual module.

**Why:** The migrator targets this external module but does not own it. The external developer has indicated willingness to accept API additions (e.g., `mf_to_double` for I/O conversion) proposed by the migrator project.

**How to apply:** When designing migrator features that require module API (new operators, conversion functions, constants), propose the API shape in the plan — the external developer will likely accept. Don't assume the module needs to be built from scratch.
