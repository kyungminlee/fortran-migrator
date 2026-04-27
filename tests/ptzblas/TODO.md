# tests/ptzblas — TODO / follow-ups

## Coverage gaps

- No tests yet for kind10 / multifloats targets — wrappers (`target_kind10`,
  `target_multifloats`) are present but have not been built. Verify by
  staging with `--target kind10` or `--target multifloats` and running
  `ctest -R '^ptzblas_'`.
- Routines present in upstream PTZBLAS but not yet exercised:
  `qmmcadd`, `qmmtcadd`, `xmmcadd`, `xmmtcadd` (mostly identical to
  the non-`c` variants in the real domain; trivial to add when needed),
  `qtzpadcpy`, `xtzpadcpy` (copy-and-pad — needs an extra B operand).
- `drshft` / `dcshft` are tested with positive offsets only. Add a
  negative-offset case once the corresponding shape conventions are
  confirmed.
- `dtzpad` / `ztzpad` exercise IOFFD = 0 only. Off-diagonal cases
  (positive/negative IOFFD, with `M ≠ N`) would be valuable for
  catching index-arithmetic bugs in the migrated code.
- Mixed-precision routine `qxvasum` is wired up via the wrapper's
  `target_dzvasum` and tested. There is no analogous `xqvasum` in the
  source set, so no symmetric test was written.

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
