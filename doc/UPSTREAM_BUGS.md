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
