# Upstream bugs in vendored Netlib sources

This document catalogues bugs found in the vendored upstream sources
(`external/lapack-3.12.1/`, `external/scalapack-2.2.3/`, etc.) that
the migrator works around without editing `external/`. Each entry
records the symptom, root cause, and the in-tree workaround. Entries
that have been reported to upstream link the tracking issue.

## How fixes are carried

Recipes accept a ``source_overrides`` field (see ``recipes/README.md``)
that maps an upstream filename to a replacement source written in
upstream shape (``DOUBLE PRECISION`` types, ``pd*``/``pz*`` symbol
names, ``dgemm`` call sites, …). The replacement goes through the
normal migration pipeline, so a single override produces correctly
renamed/promoted output for every target. The standard-precision
archive built from the unmodified ``external/`` tree is unaffected
— only the migrated extended-precision archive carries the fix.

When the convergence picker would otherwise pick the un-fixed C/S
half over the patched D/Z half, the recipe's ``prefer_source`` field
pins the correct canonical (the rank picker doesn't recognize
ScaLAPACK's ``pd*``/``pz*`` two-letter prefix — first character is
always ``P`` — so it sorts alphabetically by file name).

---

## ScaLAPACK 2.2.3: `p?lanhs.f` NPROW=1 underestimate (1/F/I norms)

**Symptom.** Migrated `pqlanhs` / `pxlanhs` (and the upstream halves
they came from) return 1-norm, Frobenius-norm, and infinity-norm
values 10–20% smaller than the reference for upper-Hessenberg matrices
of size n ≥ 32 with MB ≥ 8. The max-element norm (`'M'`) appears to
pass — but only by luck on random matrices. The error is independent
of process count: it reproduces with `mpirun -np 1` as well as 2×2 grids.

**Root cause.** The NPROW=1 first-block code fails to advance the
local row counter `II` after processing the first block of columns.
The inner-loop bound `MIN(II + LL - JJ + 1, IIA + NP - 1)` is
supposed to stop at the local row of column `LL`'s subdiagonal —
which depends on `II` tracking the local row corresponding to the
top of the current column block. The structure is:

```fortran
IF( NPROW.EQ.1 ) THEN
   IF( MYCOL.EQ.IACOL ) THEN
      DO LL = JJ, JJ+JB-1
         ...inner loop bounded by II+LL-JJ+1...
      END DO
      JJ = JJ + JB             ! JJ advances
   END IF
   IACOL = MOD( IACOL+1, NPCOL )
   ! II is *not* advanced here — bug.

   DO J = JN+1, JA+N-1, NB
      ...inner loop using stale II...
      JJ = JJ + JB
      II = II + JB             ! main loop advances II every iteration
   END DO
END IF
```

When control enters the main loop for the second column block:

* `JJ` has been advanced to `JN + 1` (correct)
* `II` is still `IIA` (should be `IIA + JB`)

So for the first column of block 2 (`LL = JN + 1`):

```
II + LL - JJ + 1 = IIA + (JN+1) - (JN+1) + 1 = IIA + 1
```

That's row 2. The correct subdiagonal for column `JN + 1` is at
row `JN + 2` (e.g., row 10 for `JN = 8`). The inner loop reads only
2 rows where it should read 10 — eight elements per column are
silently dropped. After the main loop's first iteration `II` finally
advances, but it's still `JB` short of where it should be, and that
gap propagates for the rest of the matrix.

The NPROW>1 first-block code already does the analogous
`IF MYROW.EQ.IAROW THEN II = II + JB` advance after its own
first-block code, so that path is correct.

**Why M-norm passes by luck.** `MAX(|A(i,j)|)` doesn't care if you
skip elements as long as the actual maximum sits inside the rows you
do read. For random uniform entries the max element typically lands
in the upper-left, which is always in the kept range. M-norm tests
silently agreeing with the reference is not evidence that the code is
correct — it's evidence that the test inputs aren't adversarial.

The 1/F/I-norms are sums (or sums of squares); dropped elements
directly reduce the result by the missing fraction.

**Affected files.** All four precision halves carry the identical
buggy NPROW=1 path:

* `external/scalapack-2.2.3/SRC/pdlanhs.f` (used by our migrated D-half)
* `external/scalapack-2.2.3/SRC/pzlanhs.f` (used by our migrated Z-half)
* `external/scalapack-2.2.3/SRC/pslanhs.f` (S-half — bug present but
  not exercised by our extended-precision targets, since migration
  picks D as canonical for the real family)
* `external/scalapack-2.2.3/SRC/pclanhs.f` (C-half — same status)

The single-precision halves remain buggy in our standard-precision
archive (`libblas`, `liblapack`, …) since we link those directly
from `external/`. Standard-precision callers see the upstream
behavior. Only the migrated extended-precision archives carry the
fix.

**Fix.** Add the missing `II = II + JB` after the first-block code
in each of the four norm branches (M / 1 / I / F). Mirrors the
NPROW>1 branch's existing update.

```fortran
IF( NPROW.EQ.1 ) THEN
   IF( MYCOL.EQ.IACOL ) THEN
      DO LL = JJ, JJ+JB-1
         ...
      END DO
      JJ = JJ + JB
   END IF
   II = II + JB                 ! ← add this line
   IACOL = MOD( IACOL+1, NPCOL )
   ...
END IF
```

After the fix, all four norms agree with the reference to full
target precision: ~33 digits on KIND=16, ~19 on KIND=10, ~32 on
multifloats double-double.

**Workaround in tree.**

* `recipes/scalapack/source_overrides/pdlanhs.f`
* `recipes/scalapack/source_overrides/pzlanhs.f`

Wired via `recipes/scalapack.yaml`'s `source_overrides:` map.
`PDLANHS` and `PZLANHS` are pinned in `prefer_source:` so the
patched D/Z halves win convergence over the un-fixed C/S siblings.

**Why upstream Netlib's test suite never caught it.** Their tests
generate random matrices and lean heavily on M-norm coverage; the
1/F/I-norm paths on truly upper-Hessenberg input simply aren't
exercised. The disagreement only surfaces when you generate a
genuinely Hessenberg matrix (zero below the subdiagonal) and compare
the sum-norms to a serial reference — exactly what
`tests/scalapack/auxiliary/test_p[dz]lanhs.f90` do.

**Test drivers.**

* `tests/scalapack/auxiliary/test_pdlanhs.f90` — real Hessenberg, all four norms.
* `tests/scalapack/auxiliary/test_pzlanhs.f90` — complex Hermitian-Hessenberg, all four norms.

Both PASS to full target precision on all three targets after the fix.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `p?lanhs.f` IAROW double-advance (NPROW>1)

**Symptom.** A sibling bug to the NPROW=1 one above, hiding in the
NPROW>1 branch of the same routines. With a 2×2 (or any NPROW>1)
grid, the 1/F/M norms of an upper-Hessenberg matrix come out 20–50%
under-reported; I-norm sometimes matches by chance because every row
sum is contributed exactly once (just by the wrong rank). With
NPROW=1 the bug is silent because `MOD(_, 1) == 0` collapses the
miscompute.

**Root cause.** The first-block + per-block-iteration code maintains
a tracked owner row `IAROW` that should advance one row per block.
The upstream pattern is:

```fortran
INXTROW = MOD( IAROW + 1, NPROW )    ! before the first-block code
...first-block work...
IAROW = INXTROW                      ! step IAROW forward by one block
IAROW = MOD( IAROW + 1, NPROW )      ! WRONG: advances IAROW *again*
```

`INXTROW` is computed once before the first-block code and never
recomputed inside the main loop, so the second assignment is a typo
for `INXTROW = MOD( INXTROW + 1, NPROW )`. As written, `IAROW` skips
one row owner per iteration; for NPROW=2 this leaves `IAROW` stuck
on row 0 forever, so half the block rows are read by the "this row
owns the diagonal element" branch on the wrong rank. The 1/F norms
sum the dropped contributions to zero and underestimate; the M-norm
similarly misses entries.

The pattern repeats eight times in each of `pdlanhs.f` / `pzlanhs.f`:
once per norm (M, 1, I, F) × {first-block code, main loop body}.

**Affected files.** Same set as the NPROW=1 bug above — all four
precision halves (`p[sdcz]lanhs.f`) carry the identical pattern.

**Fix.** Replace the duplicate `IAROW = MOD(IAROW+1, NPROW)` with
`INXTROW = MOD(INXTROW+1, NPROW)` so `INXTROW` advances each
iteration and IAROW correctly tracks the next-block owner:

```fortran
IAROW   = INXTROW
INXTROW = MOD( INXTROW + 1, NPROW )
```

**Workaround in tree.** Same pair as the NPROW=1 fix:

* `recipes/scalapack/source_overrides/pdlanhs.f`
* `recipes/scalapack/source_overrides/pzlanhs.f`

Both the NPROW=1 `II = II + JB` patch and this IAROW fix live in the
same override file.

**Why upstream's tests miss it.** The upstream LIN driver feeds
random general matrices and matches against `DLANGE`, not against a
hand-computed reference for genuinely upper-Hessenberg input. Since
`PDLANHS` returns a value that's a slight underestimate on random
matrices, the driver's "result is finite and not too far from
DLANGE's value" assertion accepts it. Only when the matrix is
zeroed below the first subdiagonal (so the inner-loop bound
`II + LL - JJ + 1` actually constrains anything) does the bug
surface, which is what our `test_p[dz]lanhs.f90` drivers do.

**Upstream report.** Not yet filed.

---

## ScaLAPACK 2.2.3: `p?geequ.f` column-scale reduction wrong axis

**Symptom.** `PDGEEQU`/`PZGEEQU` return correct row scaling and
`AMAX` on multi-MYCOL grids, but `COLCND` and the `INFO` zero-column
detection are computed per grid column instead of globally. For a
2×2 grid, ranks `(_,0)` see a `COLCND` derived from columns 0..NB-1,
2*NB..3*NB-1, ...; ranks `(_,1)` see `COLCND` from the other set.
With our `tests/scalapack/auxiliary/test_p[dz]geequ.f90` random
matrices the disagreement is typically ~5% on `COLCND` and on
individual `C(j)` values that depend on the global RCMAX. The
real-precision test passes by chance on the seeds we use; the
complex test fails consistently.

**Root cause.** The C scale-factor computation has two reductions:

```fortran
! Each rank computes max over its local column-block:
DO J = JJA, JJA+NQ-1
   DO I = IIA, IIA+MP-1
      C(J) = MAX( C(J), CABS1( A(IOFFA+I) ) * R(I) )
   END DO
END DO
! Combine per-column maxes across grid rows (same MYCOL):
CALL DGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, NQ, C(JJA), 1, ...)

! Now compute scalar RCMAX/RCMIN over the local NQ columns:
DO J = JJA, JJA+NQ-1
   RCMAX = MAX( RCMAX, C(J) )
   RCMIN = MIN( RCMIN, C(J) )
END DO
! Combine across processes — this is the bug:
CALL DGAMX2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, RCMAX, 1, ...)
CALL DGAMN2D( ICTXT, 'Columnwise', COLCTOP, 1, 1, RCMIN, 1, ...)
```

After the first `'Columnwise'` reduction on `C(JJA..)`, every rank in
a given grid column already holds the max-folded `C` values for its
column block. The second reduction (on the scalar `RCMAX`) is also
`'Columnwise'` — combining across `MYROW`-with-same-`MYCOL`, ranks
that already agree. So it's a no-op, and ranks holding *different*
column blocks (different `MYCOL`) keep their disjoint local
extrema, leading each grid column to compute `COLCND` from its own
block max/min instead of the global. The `IGAMX2D` for `INFO` (zero
column detection) has the same axis mistake.

The mirror axis would be `'Rowwise'` with `ROWCTOP`, which combines
across `MYCOL`-with-same-`MYROW` — exactly what's needed to merge
the disjoint per-grid-column blocks into a single global value.

The R scale-factor reductions are correct: after a `'Rowwise'` C(R)
reduce across `MYCOL`, the scalar RCMAX reduce uses `'Columnwise'` to
combine across `MYROW`, which is the right mirror.

**Affected files.**

* `external/scalapack-2.2.3/SRC/pdgeequ.f`
* `external/scalapack-2.2.3/SRC/pzgeequ.f`
* `external/scalapack-2.2.3/SRC/psgeequ.f` (same bug, untested by us)
* `external/scalapack-2.2.3/SRC/pcgeequ.f` (same bug, untested by us)

Lines (in `pdgeequ.f`): 332, 334, 346 — the three calls to
`DGAMX2D / DGAMN2D / IGAMX2D` after the C(JJA) accumulation block.
The same three lines exist in `pzgeequ.f` (with `DGAMX2D` /
`DGAMN2D` / `IGAMX2D` — note the `D`, not `Z`, for the real-typed
extrema in the complex routine).

**Fix.** Change `'Columnwise'` to `'Rowwise'` and `COLCTOP` to
`ROWCTOP` on those three calls. The patched override in
`recipes/scalapack/source_overrides/p[dz]geequ.f` carries the fix
plus an inline comment block explaining the mirror.

**Why upstream's tests miss it.** The driver compares row/column
scalings against an LU-based condition-number proxy, not against
sequential `DGEEQU` directly. The block-local `COLCND` happens to
fall in the same "no scaling needed" bucket for random matrices,
so the no-op path through the rest of the driver agrees. A
seed-by-seed numerical comparison surfaces the disagreement
immediately, which is what `test_p[dz]geequ.f90` do.

**Upstream report.** Not yet filed.

## ScaLAPACK 2.2.3: `p?posvx.f` LWMIN too small (PDPOCON / PZPOCON aborts)

**Symptom.** Calling migrated `pqposvx` / `pxposvx` (or the upstream
`PDPOSVX` / `PZPOSVX`) in single-thread or any nontrivial grid
configuration aborts inside the internal `PDPOCON` / `PZPOCON` call
with `On entry to P[QX]POCON parameter number 10 had an illegal value`.
The outer `*POSVX` workspace query returns `WORK(1) = 3*DESCA(LLD_)`,
which the caller dutifully allocates — but `PDPOCON` then enforces a
much larger `LWMIN` and aborts.

**Root cause.** `pdposvx.f:430` and `pzposvx.f:429` set
``LWMIN = 3*DESCA( LLD_ )`` where the documented contract (line 311 of
either file's prologue) is

```
LWORK = MAX( PDPOCON( LWORK ), PDPORFS( LWORK ) )
```

`PDPOCON`'s actual `LWMIN` is

```
2*NPMOD + 2*NQMOD +
MAX( 2, MAX( NB*ICEIL(NPROW-1, NPCOL),
             NQMOD + NB*ICEIL(NPCOL-1, NPROW) ) )
```

and `PDPORFS`'s is `3*NPMOD`. For any nontrivial grid (e.g. 2×2 with
`NB ≥ 8`), the documented formula exceeds `3*DESCA(LLD_)`. Same shape
for `PZPOSVX` with `LWMIN_PZPOCON = 2*NPMOD + MAX(2, …)` (no `2*NQMOD`
term — that part lives in `LRWMIN`) and `LWMIN_PZPORFS = 2*NPMOD`. The
`LRWMIN = MAX( 2*NQ, NP )` formula in `pzposvx.f:430` is also too small
versus `LRWMIN_PZPOCON = 2*NQMOD` and `LRWMIN_PZPORFS = NPMOD`.

**Files affected.**
- `external/scalapack-2.2.3/SRC/pdposvx.f` (line 430).
- `external/scalapack-2.2.3/SRC/pzposvx.f` (lines 429–430).

**Fix.** Patched overrides in
`recipes/scalapack/source_overrides/p[dz]posvx.f` recompute
`NPMOD`/`NQMOD` (PDPOCON's unadjusted NUMROCs), then set
`LWMIN = MAX( PDPOCON_LWMIN, PDPORFS_LWMIN )` (and `LRWMIN` accordingly
for the complex variant). `recipes/scalapack.yaml` declares the
overrides plus the matching `prefer_source: PDPOSVX, PZPOSVX` pin to
keep the canonical-rank picker from selecting the un-fixed S/C halves.

**Companion bug — `pzposvx.f` LRWMIN too small.** The same file's
`LRWMIN = MAX( 2*NQ, NP )` doesn't cover `PZLANHE('1', ...)`'s
documented `RWORK` requirement of `2*Nq0 + Np0 + LDW`, which `PZPOSVX`
calls at line 610 to compute ANORM. The shortfall doesn't abort
immediately — `PZLANHE` writes past the end of the caller's RWORK and
corrupts the heap, surfacing as a `free()` SIGSEGV during cleanup or a
few iterations later. Fix: `LRWMIN >= 2*NQMOD + NPMOD + LDW` (LDW=0
on square grids; we use `NPMOD` as a safe upper bound).

**Why upstream's tests miss it.** The shipped `*posvx` test drivers
allocate `WORK` / `RWORK` from precomputed bounds that exceed the
flawed `LWMIN` / `LRWMIN`, so the in-tree workspace queries are never
the binding constraint. A driver that trusts the queries (the natural
call pattern) reproduces the LWORK abort immediately and the LRWORK
heap corruption a few iterations later.

**Upstream report.** Not yet filed.

---
## ScaLAPACK 2.2.3: `pdtrsen.f` IWORK(1:N) early-write during LQUERY

**Symptom.** Calling `PDTRSEN` with `LIWORK = -1` (workspace query)
and a 1-element `IWORK` corrupts the heap, surfacing as a
`free(): invalid pointer` abort on the next `deallocate` of any
unrelated buffer. The query *also* returns the correct `LIWMIN` in
`IWORK(1)`, so a follow-up full call with a properly sized `IWORK`
runs to completion and produces a numerically correct eigenvalue
spectrum — the heap damage shows up only at cleanup.

**Root cause.** The `LQUERY` gate at `pdtrsen.f:475` enters the
parameter-validation block on either `INFO=0 OR LQUERY`. Inside that
block, the `SELECT(1:N) → IWORK(1:N)` integer-conversion loop at
lines 499-525 writes `IWORK(K)` (and `IWORK(K+1)` for 2-by-2 block
boundaries) for `K = 1..N`. Line 537 then broadcasts `IWORK(1:N)` via
`IGAMX2D(...,IWORK,N,...)`. Only after all of this is the `LQUERY`
return reached at line 619-622, which sets `IWORK(1) = LIWMIN`. The
caller passing `IWORK(1)` for the query (the documented contract) has
`IWORK(2..N)` written past the end of its 1-element buffer.

**Files affected.**
- `external/scalapack-2.2.3/SRC/pdtrsen.f` (lines 499-538).

**Fix.** Wrapper-side, mirroring the `target_pdsyevx` `WORK_T(3)` bump
for upstream `PDSYEVX`'s `WORK(1:3)` early-write. `target_pdtrsen` in
`tests/scalapack/common/target_scalapack_body.fypp` now allocates a
local `iwork_t(max(1,n))` for the LQUERY branch and forwards it to
upstream instead of the caller's `IWORK`; `iwork_t(1)` (= LIWMIN) is
copied back to `iwork(1)` after the query returns.

**Why upstream's tests miss it.** Most callers (including ScaLAPACK's
own test harness) allocate `IWORK` to a precomputed `LIWMIN` upper
bound *before* the query, so the early-write fits. A caller that
follows the documented two-pass query contract (`IWORK(1)` for the
query, then size to `LIWMIN`) trips the heap corruption.

**Upstream report.** Not yet filed.
