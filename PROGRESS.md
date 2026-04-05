# Progress / TODO

## LAPACK convergence status

LAPACK 3.12.1 type-migration currently produces **349 divergent
pairs** out of 1018 migrated routines (down one after closing the
rename-pattern cache bug below). All remaining entries are
genuine upstream source-level S/C vs D/Z drift — see `NOTE.md`
for full categorization. The D/Z version is retained on disk as
canonical for each pair; no further migrator work is expected to
reduce this count without risking hiding real semantic
differences.

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

The current 349 is down from earlier runs via these pipeline
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

## Fortran `converge` status (BLAS + LAPACK)

After the session's fixes, the authoritative post-migration
report (`converge RECIPE OUTPUT_DIR`) produces:

| library | diverged pairs | remaining category                  |
|---------|----------------|--------------------------------------|
| BLAS    | 20             | local-var / dummy-arg prefix drift   |
| LAPACK  | 464            | SROUNDUP_LWORK asymmetry, REAL(int,KIND=N) casts, ILAPREC 'D'/'E', USE ISO_FORTRAN_ENV REAL32/REAL64, algorithmic drift |

BLAS's 20 entries are entirely the S/D/C/Z prefix convention on
local variables (`SX`/`DX`, `STEMP`/`DTEMP`, `CA`/`ZA`, …). These
are not in the symbol table and are outside the migrator's reach
without a local-rename policy.

LAPACK's 464 are dominated by the documented categories (see
`NOTE.md`). The `--grep` / `--exclude` flags on `converge` are
the triage knobs.

### Fixes this session (uncommitted)

1. **Rename-pattern cache bug** (`pyengine/fortran_migrator.py`).
   `_RENAME_PATTERN_CACHE` was keyed by `id(rename_map)`. When
   `_migrate_with_flang` built a fresh `file_rename_map` per file,
   Python re-used the freed dict's address, causing cache hits
   against a stale pattern — roughly 30 subroutine/call names per
   LAPACK run kept their original prefix. Now keyed by
   `frozenset(upper_map.items())`. BLAS `migrate` divergences
   dropped from 33 → 3.
2. **Bare-literal canonicalization** (`replace_literals`). Bare
   `0.0`/`1.0` now get `E0_{kind}` appended so C-source
   `(0.0,0.0)` converges with Z-source `(0.0d0,0.0d0)`. String
   literals (FORMAT specs) and Fortran logical operators
   (`.EQ.`/`.AND.`/…) are excluded via string-segmentation and
   operator masking. `E+0` exponents are also normalized to `E0`.
3. **Light normalizer** (`_light_normalize`). Now sorts
   `INTRINSIC`/`EXTERNAL` argument lists and bare-identifier
   lists after simple type specs (`REAL(KIND=N) A,B,C`,
   `INTEGER A,B,C`), which are purely cosmetic upstream ordering
   drift. Per-line stripping is applied before the sort regex so
   the `^` anchor matches fixed-form leading whitespace.
4. **Language gate on converge**. Non-Fortran recipes now print
   a clear skip message and return an empty report rather than
   silently reporting 0 diverged.

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

## Mid priority

- **C support in `converge`**: currently only Fortran recipes are
  verified. BLACS (pure C) and PBLAS (C + Fortran kernels via
  `extra_symbol_dirs`) get a skip message. Needs a per-file
  in-memory C migration function (the existing `migrate_c_directory`
  is directory-level cloning), plus a light C normalizer
  (whitespace/comment-only) and a C pair discovery pass over
  `source_dir`. The PBLAS Fortran kernels live in subdirs
  (`PBBLAS`, `PTZBLAS`) that the current pair-finder doesn't
  walk — pair discovery should follow the same transitive
  `extra_symbol_dirs` logic used for symbol scanning.

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
