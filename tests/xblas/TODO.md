# tests/xblas — TODO

## Coverage status

| Target       | Status    | Notes                                |
|--------------|-----------|--------------------------------------|
| kind16       | 43/43 ✓  | Quad working precision               |
| kind10       | 43/43 ✓  | 80-bit extended; built+passes clean  |
| multifloats  | 43/43 ✓  | DD-on-double; 85/85 case-level PASS  |

All **43 family heads** of XBLAS are tested and passing on the kind
targets — full user-facing surface coverage.  See `level1/`, `level2/`, `level3/` for the test
drivers.

```
level1: axpby, dot, sum, waxpby                                  — 8 tests (real & complex)
level2: gemv, gemv2, gbmv, gbmv2, ge_sum_mv,
        symv, symv2, sbmv, spmv,
        hemv, hemv2, hbmv, hpmv,
        trmv, tpmv, tbsv, trsv                                  — 32 tests
level3: gemm, symm, hemm                                          — 5 tests
                                                                  ──────────
                                                                  43 total
```

The complex-symmetric variants (`zsymv_x`, `zsbmv_x`, `zspmv_x`) and
the head/tail (`*2`) variants and `ge_sum_mv_x` have no standard-BLAS
analogue; their references are hand-coded at REAL(KIND=16) precision
in `common/ref_quad_xblas.f90`.

## Open follow-ups

### Other targets — RESOLVED

All three targets (kind10 / kind16 / multifloats) build and pass
without per-target adjustment. Confirmed 2026-04-30: 43/43 tests
on multifloats, 85/85 case-level PASS. The wrapper module's
fypp parameterization on the prefix-character pair was sufficient;
no header_patch / tolerance bumps required for the dd-on-double
target.

### Mixed-input variants

The 402 mixed-input-precision variants (`BLAS_dgemv_d_s_x`,
`BLAS_zaxpby_c`, ...) are skipped by the recipe — they're not
expressible in a single extended-precision target.  If a future
caller needs one of these, the path is to add a hand-written
override in `recipes/xblas/mfc_overrides/` rather than trying to
push them through the family classifier.

### LAPACK xx-family unblocking (P17xx) — RESOLVED

`tests/lapack/TODO.md` §P17xx documented 18 LAPACK extra-precise
iterative-refinement drivers (`gesvxx`, `gbsvxx`, `gerfsx`, ...)
blocked on XBLAS being available.  Now landed — 18/18 PASS on all
three targets.  Resolution: `${LIB_PREFIX}lapack` PUBLIC-links
`${LIB_PREFIX}xblas` so migrated qla/ela/mla `*_extended` objects
resolve via qxblas/exblas/mxblas, and a quad-precision Fortran
bridge (`tests/lapack/reflapack/refxblas_quad_bridge.f90`) supplies
the upstream-named `blas_dgemv_x_` / `blas_zhemv_x_` / etc. that
reflapack_quad's `dla_*_extended.f` calls — forwarding to the
quad-promoted refblas_quad standard L2 routines.  See
`tests/lapack/TODO.md` §P17xx for the full write-up.

## Cross-tests/ items deferred

None.  `tests/CMakeLists.txt` already iterates `STAGED_LIBRARIES`,
so `tests/xblas/CMakeLists.txt` is picked up automatically once
`xblas` appears in `LIBRARY_ORDER` (which it now does).

## Known XBLAS oddities to remember

- All uniform-precision XBLAS entry points carry the `_x` suffix and
  the `prec` enum parameter — there is no `BLAS_dgemv` (without
  `_x`) in the source distribution.  Tests pass `blas_prec_extra=214`
  via the wrapper.
- The `dot_x` family computes `r := α·(x.y) + β·r` (with `r` as
  inout), not the bare `(x.y)` of standard `ddot`.
- The XBLAS conj enum is `blas_no_conj=192`, `blas_conj=191` —
  *not* the trans-style 111/112 numbers — get this wrong and only
  the conj=1 case fails (conj=0 is identical to no-op conjugation).
- The XBLAS triangular routines (`trmv_x`, `tpmv_x`, `tbsv_x`,
  `trsv_x`) take an `α` argument that standard BLAS triangular
  routines do not.  Reference implementations call standard
  `*trmv`/`*trsv` then scale by α (or scale before, for solves).
- Complex parameters in the C API are `void *`; from Fortran via
  ISO_C_BINDING we declare them as `complex(ep)` (passed by
  reference).  Memory layout matches because `complex(KIND=16)` is
  two adjacent `__float128` values.
- The XBLAS `f2c-bridge` macros default to no-mangling; the migrator
  build sets `-DCONFIG_FC_UNDERSCORE` so the Fortran-callable names
  end with `_` (matching gfortran's emit convention).  If we ever
  add an Intel / Cray Fortran build we'll need a different gate
  (`-DCONFIG_FC_UCASE` or `-DCONFIG_FC_DBL_UNDERSCORE`).
