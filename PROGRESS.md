# Progress / TODO

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
