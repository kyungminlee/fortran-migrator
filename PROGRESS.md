# Progress / TODO

## LAPACK convergence status

LAPACK 3.12.1 type-migration currently produces **350 divergent
pairs** out of 1018 migrated routines. All are genuine upstream
source-level S/C vs D/Z drift — see `NOTE.md` for full
categorization. The D/Z version is retained on disk as canonical
for each pair; no further migrator work is expected to reduce
this count without risking hiding real semantic differences.

Breakdown:

| category                    | count |
|-----------------------------|-------|
| `SROUNDUP_LWORK` asymmetry  | 237   |
| `ILAPREC` `'D'` vs `'E'`    | 9     |
| label renumbering           | ~20   |
| unused `PARAMETER(ONE=…)`   | ~15   |
| expression formatting drift | ~15   |
| `INTRINSIC` list mismatch   | 3     |
| algorithmic drift           | ~5    |
| workspace-query dummy args  | 1     |
| misc                        | ~45   |

### Recent canonicalization improvements

The current 350 is down from earlier runs via these pipeline
changes (all in `pyengine/pipeline.py`):

- `_strip_real_cmplx_casts()` — strips `REAL(expr, KIND=N)`,
  `CMPLX(expr, KIND=N)`, 3-arg `CMPLX(re, im, KIND=N)`, and
  single-arg `REAL(X)` (but not single-arg `CMPLX(X)` which
  changes type).
- `_has_top_level_operator()` — only wraps stripped expressions
  in parens when they contain a top-level `+`/`-`, since `*`/`/`
  are associative across adjacent operators.
- Numeric literal collapse: `1.0` / `1.` / `1` → `1`.
- `(KIND=N)` suffix stripping on `REAL`/`COMPLEX` type-specs
  (ordered **before** the `@`-collapse so it still matches).
- `iso_fortran_env` kind unification: `real32`/`real64`/`real128`
  → `realK`.
- Case-insensitive uppercase; `END IF` → `ENDIF` collapse;
  `IMPLICIT NONE` stripping.
- Classifier `_size` tiebreaker: added `min(positions)` key to
  prevent DSYTRS/DSYTRD-style false pairings (audit: 0
  suspicious pairings in LAPACK).

## BLAS convergence status

3 divergent pairs (all documented in `NOTE.md`): `sdot/ddot`,
`scasum/dzasum`, `srotmg/drotmg`. All benign.

## Low priority

- **Two-step classifier (quartet-first)**: As a defensive layer on top
  of the existing pair-finder, do a first pass that treats S/D/C/Z as
  interchangeable (replace with `#T`) and greedily assign fully-populated
  quartets before running the current pair-discovery algorithm on the
  remainder. Not needed for any currently-known case — the
  min-position tiebreaker already handles the DSYTRS/DSYTRD-style
  false pairings, and a full audit of LAPACK found zero suspicious
  pairings — but would provide extra robustness against future
  libraries with quirky prefix conventions.
