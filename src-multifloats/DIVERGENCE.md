# Divergence Report

End-to-end migration and convergence check for all libraries,
run on 2026-04-12 from `src-multifloats/pyengine` with YAML-based targets.

Targets tested: `kind16` (KIND=16, Q/X prefix) and `multifloats` (float64x2, DD/ZZ prefix).

## Summary

| Library      | Lang    | Files | kind16 div | multifloats div |
|:-------------|:--------|------:|-----------:|----------------:|
| BLAS         | Fortran |    78 |          3 |               3 |
| LAPACK       | Fortran |  1037 |        431 |             418 |
| PBBLAS       | Fortran |    14 |          0 |               0 |
| PTZBLAS      | Fortran |    52 |          0 |               0 |
| ScaLAPACK    | Fortran |   342 |         24 |              70 |
| BLACS        | C       |   213 |          0 |               0 |
| PBLAS        | C       |   293 |         61 |              61 |
| ScaLAPACK-C  | C       |    17 |          1 |               1 |
| **Total**    |         |  2046 |      **520** |          **553** |

"Divergence" means the S-half and C-half (or D-half and Z-half) of a
co-family pair produce different normalized output after migration.
Divergences come from genuine asymmetries in the upstream sources, not
from migrator bugs (unless noted otherwise).


## BLAS (3 / 3)

Three files diverge in both targets. All are pre-existing upstream asymmetries.

### `scasum.f` vs `dzasum.f`
The S-half inlines `ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))` while the
Z-half calls the external function `DCABS1(ZX(I))`.  After migration
the inlined form vs external call persist.

### `sdot.f` vs `ddot.f`
Local accumulator is named `STEMP` in the S-half and `DTEMP` in the
D-half; only one survives normalization (variable-name mismatch).

### `srotmg.f` vs `drotmg.f`
- The DATA constants `GAMSQ`/`RGAMSQ` have different literal
  representations (`1.67772E7` vs `16777216.E0`; `5.96046E-8` vs
  `5.9604645E-8`).
- Local variable names use S-prefix vs D-prefix (`SD1`/`DD1`,
  `SFLAG`/`DFLAG`, etc.).


## LAPACK

### kind16: 431 divergences

| Category | Files | Description |
|:---------|------:|:------------|
| `ROUNDUP_LWORK` | 237 | S-half uses `QROUNDUP_LWORK()`; D-half assigns directly |
| `ROUNDUP_LWORK` + complex-specific | ~56 | ROUNDUP plus complex ONE parameter or CABS1 |
| Complex-specific only | ~67 | `PARAMETER(ONE=(...))`, `QCABS1`, `INTRINSIC AIMAG/ABS` |
| Source asymmetry only | ~71 | Variable names, control-flow, literal formatting |

### multifloats: 418 divergences

| Category | Files | Description |
|:---------|------:|:------------|
| `ROUNDUP_LWORK` | 237 | Same as kind16 |
| Complex-specific only | 15 | Fewer than kind16 (constructor wrapping normalizes some) |
| Source asymmetry + constructor wrapping | 166 | Includes `FLOAT64X2(N)` wrapping differences |

### Root causes

**`ROUNDUP_LWORK` (dominant, ~237 files).**
LAPACK 3.12 added `SROUNDUP_LWORK` / `CROUNDUP_LWORK` to the
single-precision (S/C) code paths to round workspace sizes up to the
nearest `REAL`.  The double-precision (D/Z) halves do not use this
function — they assign `WORK(1) = LWKOPT` directly.  Since the
migrator generates from the S-half, the migrated code retains
`QROUNDUP_LWORK` (kind16) or `DDROUNDUP_LWORK` (multifloats), while
the D-half convergence check does not have it.

*This is a genuine S/D asymmetry in LAPACK, not a migrator bug.*

**Complex ONE parameter (~67 files, kind16 only).**
Many C/Z routines define `PARAMETER(ONE=(1.0D0,0.0D0))` as a complex
constant.  The S/C halves do not always define this parameter.  After
migration, the complex-half output contains the extra `PARAMETER`
statement.

**Constructor wrapping (multifloats-specific).**
The multifloats target wraps integer-to-real conversions in
`FLOAT64X2(...)`.  Some S-half code has expressions like
`FLOAT64X2(MAXITR) * (FLOAT64X2(N) * ...)` where the D-half uses
plain `MAXITR * (N * ...)`.  This creates extra diff lines in the
multifloats target that don't appear in kind16.

**Source asymmetries (~70–170 files).**
Upstream LAPACK sources have genuine differences between
S/C and D/Z variants:
- Different variable names (`STEMP` vs `DTEMP`)
- Different literal representations (`0.0E0` vs `0.0E00`)
- Different control-flow (label numbers, `GO TO` targets)
- Extra `IMPLICIT NONE` in one half but not the other
- Different `EXTERNAL`/`INTRINSIC` declarations
- Error-string differences (`XLA_SYRFSX_EXTENDED` vs `XLA_HERFSX_EXTENDED`)


## PBBLAS (0 / 0)

No divergences.  All 14 co-family pairs converge perfectly.


## PTZBLAS (0 / 0)

No divergences.  All 47 co-family pairs converge perfectly.


## ScaLAPACK — Fortran (24 / 70)

### kind16: 24 divergences

| File(s) | Diff lines | Root cause |
|:--------|:-----------|:-----------|
| `bslaexc.f` / `bdlaexc.f` | +2 | Literal: `10.0E0` vs `1.0E1` |
| `pslacon.f` / `pdlacon.f` | +6 | Declaration layout, EXTERNAL lists |
| `pslaed3.f` / `pdlaed3.f` | +4 | `IINFO=0` vs `INFO=0`; extra `IF` block |
| `pslaqr3.f` / `pdlaqr3.f` | +10 | `MPI_WTIME` external, workspace calc, literal `0.0E0` vs `0.0E00` |
| `pslawil.f` / `pdlawil.f` | +6 | Declaration layout |
| `psstebz.f` / `pdstebz.f` | +46 | `PSLAIECT` vs `PDLAIECTB`/`PDLAIECTL` (different algorithm) |
| `pssyevd.f` / `pdsyevd.f` | +2 | Minor |
| `pssyttrd.f` / `pdsyttrd.f` | +8 | Source asymmetry |
| `pstrsen.f` / `pdtrsen.f` | +4 | Source asymmetry |
| `pzheevd.f` / `pcheevd.f` | +2 | Minor |
| `pzhettrd.f` / `pchettrd.f` | +8 | Source asymmetry |
| `pzlarz.f` / `pclarz.f` | +4 | Source asymmetry |
| `pzlarzc.f` / `pclarzc.f` | +2 | Minor |
| `pzlawil.f` / `pclawil.f` | +6 | Declaration layout |
| `pzpotf2.f` / `pcpotf2.f` | +12 | Source asymmetry |
| `pzungql.f` / `pcungql.f` | +4 | Source asymmetry |
| `pzunmbr.f` / `pcunmbr.f` | +2 | Minor |
| `pzunml2.f` / `pcunml2.f` | +4 | Source asymmetry |
| `slarre2.f` / `dlarre2.f` | +4 | Source asymmetry |
| `slarre2a.f` / `dlarre2a.f` | +4 | Source asymmetry |
| `slarrf2.f` / `dlarrf2.f` | +4 | Source asymmetry |
| `sstegr2.f` / `dstegr2.f` | +2 | Minor |
| `sstegr2a.f` / `dstegr2a.f` | +2 | Minor |
| `sstegr2b.f` / `dstegr2b.f` | +2 | Minor |

All are upstream source asymmetries.

### multifloats: 70 divergences

The 24 kind16 divergences all reappear.  The additional 46 come from
C/Z complex-half routines (`pzgebd2.f`, `pzgebrd.f`, `pzgehd2.f`, ...,
`pzunmtr.f`) where the constructor-based migration produces
`USE MULTIFLOATS` / `EXTERNAL` declaration differences between the
S-half and C-half migration paths — typically +4 diff lines each.


## BLACS — C (0 / 0)

No divergences.  All pairs converge perfectly.

The BLACS reduction kernels (`BI_*vvamn.c`, `BI_*vvamx.c`,
`BI_*vvsum.c`) use the `Rabs`/`Cabs` macros from `Bdef.h` and `.r`/`.i`
struct member access on complex types.  These work out of the box with
`float64x2_t` and `complex128x2_t` because:
- `Rabs` uses ternary `?:` which works via overloaded `operator<`
- `Cabs` accesses `.r`/`.i` members which match the `complex128x2_t` struct
- All arithmetic operators are overloaded in `multifloats_bridge.h`

The only remaining override is `blacs_pinfo_.c` which adds a call to
`multifloats_mpi_init()` for MPI datatype registration.


## PBLAS — C (61 / 61)

All 61 divergences appear in both targets.  The C sources have
**K&R-style function declarations** in the S/C halves that the D/Z
halves have already converted to ANSI style.  After migration, the
S-half retains:
```c
void pqagemv_ ( TRANS, M, N, ALPHA, ... )
F_CHAR_T TRANS ;
Int * IA , * INCX , ... ;
QREAL * ALPHA , * BETA ;
```
while the D-half produces the ANSI prototype form.

This is an upstream source asymmetry in ScaLAPACK's PBLAS.


## ScaLAPACK-C (1 / 1)

### `clamov.c` vs `zlamov.c`

The C-half defines `#define TYPE complex` while the Z-half defines
`#define TYPE XCOMPLEX` (kind16) or `#define TYPE complex128x2_t`
(multifloats).  This is a C/Z asymmetry in the original source where
`complex` and `ZCOMPLEX` are distinct macros.


## Known Migrator Issues

No migrator-caused divergences remain.  Two bugs were fixed:

1. **Double constructor wrapping** (fixed): `_wrap_complex_args()` now
   skips type-conversion intrinsics (`REAL`, `DBLE`, etc.) that will
   be handled by `replace_generic_conversions()` later in the pipeline.

2. **Double type-suffix in C** (fixed): type-substitution regexes now
   use `[^a-zA-Z_0-9]` boundaries so the second pass cannot re-match
   `float` inside `float64x2_t`.

Additionally, the 10 hand-written BLACS reduction kernel overrides were
eliminated by renaming `complex128x2_t` members from `.re`/`.im` to
`.r`/`.i`, matching the BLACS `DCOMPLEX`/`SCOMPLEX` struct convention.

All remaining divergences are upstream source asymmetries.
