# Notes

## Remaining BLAS convergence divergences

After type-migration, three BLAS routine pairs disagree between the
single- and double-precision halves. These are not migrator bugs —
they reflect genuine source-level differences between the S/C and D/Z
files in the upstream LAPACK 3.12.1 reference BLAS. The migrator keeps
the D/Z version on disk (our canonical output) and reports the
mismatch.

### 1. `sdot.f` vs `ddot.f` → `qdot.f`

Trivial statement-order difference in the initializer.

```fortran
! sdot.f
STEMP = 0.0e0
SDOT  = 0.0e0

! ddot.f
DDOT  = 0.0d0
DTEMP = 0.0d0
```

The two statements are reordered between halves. Semantically
identical; only statement order differs.

### 2. `scasum.f` vs `dzasum.f` → `qxasum.f`

The complex-absolute-value computation is inlined in the single
version but delegated to the helper `DCABS1` in the double version.

```fortran
! scasum.f (inline)
STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))

! dzasum.f (helper call)
STEMP = STEMP + DCABS1(ZX(I))
```

`DCABS1(z) = |Re z| + |Im z|`, so the two formulations compute the
same quantity. Upstream has simply never been harmonized; the single
side was written before a matching `SCABS1` helper was considered.

### 3. `srotmg.f` vs `drotmg.f` → `qrotmg.f`

Two differences:

**(a) Declaration order.** The `DATA`-initialized constants are
grouped differently in the local variable list:

```fortran
! srotmg.f
REAL GAM,GAMSQ,ONE,RGAMSQ,SFLAG,SH11,SH12,SH21,SH22,SP1,
     SP2,SQ1,SQ2,STEMP,SU,TWO,ZERO

! drotmg.f
DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,
                 DTEMP,DU,GAM,GAMSQ,RGAMSQ,TWO,ZERO,ONE
```

Declaration ordering has no semantic effect.

**(b) Constant literal precision.** The two halves write the shared
overflow/underflow constants at different precisions:

| constant  | srotmg.f       | drotmg.f        | value   |
|-----------|----------------|-----------------|---------|
| `GAM`     | `4096.E0`      | `4096.D0`       | 2^12    |
| `GAMSQ`   | `1.67772E7`    | `16777216.D0`   | 2^24    |
| `RGAMSQ`  | `5.96046E-8`   | `5.9604645D-8`  | 2^−24   |

`GAMSQ` and `RGAMSQ` differ beyond the 7-digit single-precision
mantissa: `16777216 = 2^24` exactly, whereas `1.67772E7` drops the
trailing `16`; similarly `5.9604645E-8 = 2^-24` while `5.96046E-8`
truncates earlier. After migration to `KIND=16`, the literals in the
S source still carry only single-precision accuracy, so the migrated
`QROTMG` produced from `srotmg.f` has slightly coarser constants than
the one produced from `drotmg.f`. The D-sourced version is retained.
