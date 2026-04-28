# tests/ptzblas — TODO / follow-ups

## Coverage gaps

- Mixed-precision routine `qxvasum` is wired up via the wrapper's
  `target_dzvasum` and tested. There is no analogous `xqvasum` in the
  source set, so no symmetric test was written.

## Multifloats `dzvasum` — fixed (root cause: complex-prefix collision)

The SEGV was a link-time symbol collision driven by multifloats's
choice of `V` as the complex prefix. PTZBLAS uses literal `V` as a
"vector" indicator in routine names (`dvasum`, `dvvdot`, `dzvasum`,
…) and BLAS uses `dz`/`sc` as the mixed real-from-complex prefix.
With C: V the slot-driven mechanical rename produced

  D[V]ASUM (real, slot R + literal V)  → T + V + ASUM = TVASUM
  D[Z]ASUM (mixed, slot R + slot C)    → T + V + ASUM = TVASUM

— two distinct upstream families both rename to TVASUM, and the
linker picks one nondeterministically. Calling the wrong TVASUM
(real-only `SUBROUTINE TVASUM(N, ASUM, X, INCX)` instead of mixed
`FUNCTION TVASUM(N, ZX, INCX)`) trampled stack memory and crashed.

**Fix**: switched the multifloats target's complex prefix from V to
W (W is unused as a literal anywhere in the BLAS / LAPACK /
ScaLAPACK / MUMPS source corpus). After re-staging multifloats,
all migrated complex-precision symbols use the `w*` prefix and
TVASUM resolves uniquely to the PTZBLAS real-only routine.

The migrator's `prefix_classifier.build_rename_map` now also raises
`RuntimeError` on any cross-family many-to-one rename collision, so
a bad prefix choice is surfaced at staging time rather than as a
runtime SEGV.

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
