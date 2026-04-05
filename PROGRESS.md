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

## ScaLAPACK convergence (in progress)

New `pyengine converge RECIPE OUTPUT_DIR` subcommand: the
authoritative post-migration check. Reads each pair's canonical
from disk, re-migrates the S/C sibling in memory, and compares
with a *light* normalizer (whitespace/case/`END KEYWORD` only). No
prefix `@`-collapse, no literal rewrites, no declaration sorting
— any remaining entry is either migrator slop or genuine drift.

The older `diverge` subcommand (heavy canonicalizer, both halves
migrated in memory) remains for preliminary exploration.

### Pipeline plumbing work done in this session (uncommitted)

- Extracted `_collect_all_symbols()` and `_scan_extra_dirs()` in
  `pipeline.py`. All three callers (`run_divergence_report`,
  `run_convergence_report`, `run_migration`) now walk transitive
  dependencies AND each dep's `extra_symbol_dirs`. Previously
  `run_migration` scanned only its own `extra_symbol_dirs`, so
  e.g. `DLAMCH`/`SLAMCH` (defined in LAPACK's `INSTALL/`) were
  missing from ScaLAPACK's rename map. Confirmed fix: `qlarre2.f`
  now has `QLAMCH` instead of a mix of `SLAMCH`/`DLAMCH`.
- Added new `run_convergence_report()` and `_light_normalize()`
  in `pipeline.py`.
- Added `converge` subcommand to `__main__.py`.
- Reverted the earlier `--from-disk` flag on `diverge` (subsumed
  by `converge`).

### Recipe fix in flight (uncommitted)

`recipes/scalapack.yaml` — added two deps that were missing:

```yaml
depends:
  - lapack.yaml
  - blacs.yaml    # SGEBR2D/DGEBR2D, SGESD2D/DGESD2D, STRRV2D/DTRRV2D, ...
  - pblas.yaml    # ScaLAPACK routines call PBLAS entry points directly
```

Without these, `converge` reports 263 divergences, most being
unrenamed BLACS communication routines (`SGEBR2D`, `SGESD2D`,
`STRSD2D`, `SGAMX2D`, `SGSUM2D`, etc.).

### Where to pick up

1. Rerun `migrate recipes/scalapack.yaml /tmp/sp-full --kind 16`
   with the updated recipe to regenerate canonicals.
2. Rerun `converge recipes/scalapack.yaml /tmp/sp-full --kind 16`
   and triage remaining categories. Expected leftover:
   - `TEN=10.0` vs `TEN=10` literal formatting
   - `INTRINSIC MAX,MIN,MOD,REAL` list ordering
   - `EXTERNAL` list ordering (declaration-reordering drift)
   - Parameter `_16` kind suffix asymmetry
   - Genuine algorithmic drift already in `diverge` report
     (PERT=4/8, MINRGP=3.0E-3/1.0E-3, TWO/FOUR*EPS, PSLAIECT
     branching, etc.)
3. Decide which of the above to handle in the migrator vs accept
   as genuine drift vs add to the light normalizer. The user's
   stated policy is: typedefs/renames are the migrator's job,
   declaration ordering and literal formatting are arguably the
   normalizer's job.
4. Commit the pipeline refactor + recipe deps once converge is
   stable.

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
