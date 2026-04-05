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


## LAPACK convergence divergences

LAPACK has **370 divergent pairs** out of 1018 migrated routines. The
bulk is explained by a handful of deliberate upstream asymmetries.

### 1. `SROUNDUP_LWORK` (237 files)

`SROUNDUP_LWORK(N)` (in `INSTALL/sroundup_lwork.f`) converts an
INTEGER workspace size to a REAL such that `INT(result) >= N`. It
exists because LAPACK returns the optimal workspace size in
`WORK(1)`, which has the same floating-point type as the matrix
data. A 32-bit REAL cannot represent every INTEGER exactly:

```
REAL(16777217)           = 16777216.0    ! truncates down — unsafe
REAL(16777217)*(1+eps)   = 16777218.0    ! rounds up — safe
```

If a caller casts `WORK(1)` back to INTEGER and allocates that many
elements, a truncated-down value yields a buffer too small.
`SROUNDUP_LWORK` multiplies by `1+eps` to force rounding up.

**Asymmetry.** The S- and C-precision routines call
`WORK(1) = SROUNDUP_LWORK(LWKOPT)`. The D- and Z-precision routines
simply write `WORK(1) = LWKOPT`. Reason: 64-bit REAL has a 53-bit
mantissa, so every INTEGER up to `2^53 ≈ 9×10^15` round-trips through
REAL exactly. No workspace will ever approach that bound, so the
double-precision side has no bug to fix and skips the helper call
entirely. Only `dgesdd.f` and `zgesdd.f` use `DROUNDUP_LWORK`
(defined alongside, and equivalent under D precision).

This is a deliberate choice by upstream LAPACK 3.12: the release
notes describe `SROUNDUP_LWORK` as fixing "a subtle bug with
returning LWORK as a Float" that only manifests in single precision,
and the D/Z halves were intentionally left bare. After migration to
KIND=16 the migrated `QROUNDUP_LWORK` is harmless (128-bit REAL has
a 113-bit mantissa, so the round-up branch never triggers), but the
two halves remain textually distinct and are flagged as divergences.
The D/Z version is retained.

### 2. `ILAPREC` precision-character argument (9 files)

(investigation pending)

### 3. `REAL32` / `REAL64` `iso_fortran_env` kinds (4 files)

(investigation pending)

### 4. Other source-level drift (~106 files)

(investigation pending)

#### Case study: `clahef.f` vs `zlahef.f` — loop-bound off-by-one

In the Bunch-Kaufman pivot-unwinding loops (both in the upper-
triangular and lower-triangular branches), the loop-continuation
check differs between the two halves:

| Source          | Upper branch (label 60) | Lower branch (label 120) |
|-----------------|-------------------------|--------------------------|
| `clahef.f` (C)  | `IF( J.LE.N ) GO TO 60` | `IF( J.GE.1 ) GO TO 120` |
| `zlahef.f` (Z)  | `IF( J.LT.N ) GO TO 60` | `IF( J.GT.1 ) GO TO 120` |

The inner `[CZ]SWAP` call is guarded by `J.LE.N` in **both** halves,
and its body starts with an unconditional `J = J + 1`. When the
C-side does its extra `J == N` iteration, the inner increment makes
`J == N+1`, the swap guard becomes false, and the iteration is a
no-op. The same holds symmetrically for the `J == 1` case in the
lower branch. So the behavioural difference is a benign off-by-one:
clahef executes one additional empty loop turn that zlahef skips.

This is genuine upstream drift — the single-complex and
double-complex halves of `LAHEF` diverged at some point and were
never resynced. The migrator cannot (and should not) paper over it;
the Z-sourced canonical is kept on disk.
