# tests/ptzblas — TODO / follow-ups

## Coverage gaps

- No tests yet for kind10 / multifloats targets — wrappers (`target_kind10`,
  `target_multifloats`) are present but have not been built. Verify by
  staging with `--target kind10` or `--target multifloats` and running
  `ctest -R '^ptzblas_'`.
- Routines present in upstream PTZBLAS but not yet exercised:
  `qmmtcadd`, `xmmtcadd` (transposed conjugate add — same shape as
  `mmtadd`; not yet wired up in the wrapper template).
- Mixed-precision routine `qxvasum` is wired up via the wrapper's
  `target_dzvasum` and tested. There is no analogous `xqvasum` in the
  source set, so no symmetric test was written.

## Multifloats `dzvasum` SEGFAULT

- **Symptom**: `test_dzvasum` segfaults on the multifloats target.
  All 27 other ptzblas tests pass on multifloats. Fixing the
  accumulator-wrapper init for the vvdot pair resolved the kind10
  NaN issue but not the multifloats SEGV here.
- **Diagnosis**: The migrated `tvvasum` calls `TVASUM(N, X, INCX)`
  which assigns through to `ASUM`. Likely a bug in the migrated
  TVASUM body or in how `q2t_c(complex(16) → cmplx64x2)` constructs
  the input vector — either way the crash is downstream of the
  wrapper, in the migrated kernel.
- **Action**: Defer; gate by skipping `test_dzvasum` on multifloats
  until the migrated tvvasum / q2t_c path is debugged.

## Potential migrator follow-ups (NOT in this subtree)

None observed during the kind16 build — every test passed on the
first iteration after the upstream-semantics fixes documented in
`common/ref_quad_aux.f90` were added (CABS1 vs Euclidean magnitude
for the complex abs-matvec family; non-circular row/column shift
semantics for `q[r|c]shft`; A-is-M×N / B-is-N×M shape for `qmmtadd`).
If ctest later reports a bug that looks like a migrator artifact (e.g.
`cmplx16` aliasing analogous to the pzher2k issue noted in
`tests/pblas/TODO.md`), record it here with a reference test program
left in place and an explanatory header comment.
