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

LAPACK has **350 divergent pairs** out of 1018 migrated routines,
breaking down as:

| category                | count |
|-------------------------|-------|
| `SROUNDUP_LWORK`        | 237   |
| `ILAPREC` character arg | 9     |
| other upstream drift    | 104   |

All 350 are genuine upstream source-level differences between the
S/C and D/Z halves of LAPACK 3.12.1 — not migrator bugs. The D/Z
version is retained on disk as canonical for each pair.

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

`ILAPREC(c)` (`SRC/ilaprec.f`) is a character→integer dispatcher
that maps a single-letter selector to one of the BLAST XBLAS
intermediate-precision enum constants:

```
'S' → BLAS_PREC_SINGLE     = 211
'D' → BLAS_PREC_DOUBLE     = 212
'I' → BLAS_PREC_INDIGENOUS = 213
'X'/'E' → BLAS_PREC_EXTRA  = 214
```

The value is passed to the XBLAS `BLAS_xGEMV_X` / `BLAS_xGEMV2_X`
calls inside the extended iterative-refinement drivers
(`[sdcz]la_*rfsx_extended.f`) to choose the residual precision.

**Asymmetry.** The single/complex drivers select `'D'` (one rung
above their base type):

```fortran
! sgerfsx.f  (and sgbrfsx / sporfsx / ssyrfsx,
!             cgerfsx / cgbrfsx / cporfsx / csyrfsx / cherfsx)
PREC_TYPE = ILAPREC( 'D' )         ! BLAS_PREC_DOUBLE = 212
```

The double/complex-double drivers select `'E'`:

```fortran
! dgerfsx.f  (and dgbrfsx / dporfsx / dsyrfsx,
!             zgerfsx / zgbrfsx / zporfsx / zsyrfsx / zherfsx)
PREC_TYPE = ILAPREC( 'E' )         ! BLAS_PREC_EXTRA = 214
```

Upstream picks whichever level is one rung above the base precision
of the driver — single uses double, double uses extra. The two
halves therefore request different XBLAS residual-precision modes
and the character literals (`'D'` vs `'E'`) survive into the
migrated output unchanged.

Nine diverging pairs: `sgbrfsx/dgbrfsx`, `sgerfsx/dgerfsx`,
`sporfsx/dporfsx`, `ssyrfsx/dsyrfsx`, `cgbrfsx/zgbrfsx`,
`cgerfsx/zgerfsx`, `cporfsx/zporfsx`, `csyrfsx/zsyrfsx`,
`cherfsx/zherfsx`.

After migration to `KIND=16` both halves are physically quad
precision; the character literal divergence is genuine upstream
source drift that the migrator cannot (and should not) paper over.
The D/Z version is retained on disk.

### 3. Other source-level drift (104 files)

These are miscellaneous places where the S/C and D/Z halves of a
routine were maintained independently and drifted apart in
token-level but semantically-benign ways. They group into roughly
six patterns:

#### 3a. DO-loop / GO TO label renumbering

The loop-label values were chosen differently between halves even
though the control structure is identical:

```fortran
! strsyl.f                  ! dtrsyl.f
DO 70 L = 1, N              DO 60 L = 1, N
  IF( L.LT.LNEXT ) GO TO 70   IF( L.LT.LNEXT ) GO TO 60
```

Largest contributors by line-count:
`strsyl/dtrsyl` (+172), `clatbs/zlatbs` (+58), `cptts2/zptts2`
(+54), `clatrs/zlatrs` (+52), `clatps/zlatps` (+48), `slatrs/dlatrs`
(+32), `slatbs/dlatbs` (+30), `slatps/dlatps` (+28),
`clapmt/zlapmt` and `slapmt/dlapmt` (+22 each), `cpttrf/zpttrf`
(+16), plus smaller label shuffles in `clahef`, `claqp2`, `clauu2`,
`claqp2rk`, `clarfx`, `cupmtr`, `slaexc`.

#### 3b. Unused `PARAMETER(ONE=…)` declaration on one side only

One half declares a constant that the code path never references:

```fortran
! cgehd2.f declares ONE; zgehd2.f does not (or vice versa):
COMPLEX ONE
PARAMETER(ONE=(1,0))
```

Pattern: `cgehd2/cgelq2/cgeql2/cgeqr2/cgeqr2p/cgerq2`,
`cunm2l/cunm2r/cunml2/cunmr2`, `sorm2l/sorm2r/sorml2/sormr2`,
`sorg2l`, `cungr2`, etc. `clahef_rk/zlahef_rk` splits a combined
`PARAMETER(ONE=…,ZERO=…)` into two single-parameter declarations.

#### 3c. `INTRINSIC` list mismatch

One half lists more intrinsics than the other, even though only the
common subset is actually called:

```fortran
! cggbak.f      INTRINSIC MAX
! zggbak.f      INTRINSIC INT, MAX

! slansf.f      INTRINSIC SQRT, ABS
! dlansf.f      INTRINSIC SQRT, ABS, MAX
```

#### 3d. Formatting / expression-level differences

Semantically identical expressions written slightly differently:

```fortran
! cheequb:  C0 = (-(T*DI)*DI+2*WORK(I)*DI-N*AVG)   ! outer parens
! zheequb:  C0 = -(T*DI)*DI+2*WORK(I)*DI-N*AVG

! ctrsyl3:  DWORK(2,1) = (2*NBB+NBA)
! ztrsyl3:  DWORK(2,1) = 2*NBB+NBA

! clarf1f:  C(2,1)
! zlarf1f:  C(1+1,1)                               ! unfolded constant

! chgeqz:   WORK(1) = 1
! zhgeqz:   WORK(1) = DCMPLX(1)                    ! explicit cast
```

#### 3e. Algorithmic drift — different code paths

A handful of pairs use genuinely different implementations:

```fortran
! slasq2:  IEEE = .FALSE.                          ! hardcoded
! dlasq2:  IEEE = (ILAENV(10,'DLASQ2','N',1,2,3,4).EQ.1)

! clarf1f calls CGER   (single-complex GER, no conjugate)
! zlarf1f calls ZGERC  (conjugate GER) — different BLAS kernel

! clahef_rk uses CGEMM
! zlahef_rk uses ZGEMMTR                           ! different routine
```

The migrator cannot (and should not) close these — they look
symmetric under our character-based canonicalizer only by accident,
but they do call different BLAS kernels.

#### 3f. Workspace-query dummy-argument asymmetry

`sorcsd/dorcsd` (+13): one half allocates a tiny
`REAL DUMMY(1)` to pass to lazy workspace-query calls; the other
passes the real output arrays (which the query-mode call never
touches). Both are correct.

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
