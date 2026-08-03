"""Co-family divergence reporting.

Migrates each source file of a precision family in memory and compares
the results: two sources that are supposed to converge onto the same
target file must produce the same text. Anything left over is a
divergence the recipe author has to look at.

This module is the join point between the two halves of the comparison
machinery — the migration pipeline (:mod:`pipeline`) and the pure text
canonicalization (:mod:`divergence`) — so it imports both at top level
and neither imports it. ``__main__`` is its only caller.
"""

from pathlib import Path

from .divergence import (
    _canonicalize_for_compare,
    _filter_expected_divergences,
    _filter_precision_drift,
    _strip_fortran_comments,
)
from .fortran_migrator import target_filename
from .pipeline import (
    build_module_rename_pairs,
    canonical_rank,
    classify_and_rename,
    migrate_parallel,
    postprocess_migrated,
)
from .prepare import prepare_recipe


def run_divergence_report(recipe_path: Path, target_mode=None,
                          project_root: Path | None = None,
                          parser: str | None = None,
                          parser_cmd: str | None = None,
                          apply_whitelist: bool = True) -> list[dict]:
    """Migrate every co-family source pair in-memory and return the
    normalized diff for each pair whose members disagree.

    Returns a list of ``{'target', 'canonical', 'other', 'diff'}``
    dicts, sorted by target filename. ``diff`` is the list of +/-
    lines (without context) from the unified diff of the two
    canonicalized texts.
    """
    config = prepare_recipe(recipe_path, project_root)

    _symbols, classification, rename_map = classify_and_rename(
        config, target_mode, project_root)

    # Group eligible source files by their target output name.
    src_files = sorted(
        p for p in config.source_dir.iterdir()
        if p.suffix.lower() in config.extensions
    )
    by_target: dict[str, list[Path]] = {}
    for p in src_files:
        stem_u = p.stem.upper()
        if stem_u in config.skip_files or stem_u in config.copy_files:
            continue
        if stem_u in classification.independent:
            continue
        by_target.setdefault(target_filename(p.name, rename_map, target_mode), []).append(p)

    # Migrate every member of every multi-member group in parallel.
    pairs: list[tuple[Path, Path]] = []
    for members in by_target.values():
        if len(members) < 2:
            continue
        members.sort(key=lambda p: (
            canonical_rank(p.stem, config.prefer_source),
            p.name,
        ))
        canonical = members[0]
        for other in members[1:]:
            pairs.append((canonical, other))

    all_paths = {p for pair in pairs for p in pair}
    module_rename_pairs = build_module_rename_pairs(config)
    results = migrate_parallel(all_paths, rename_map, target_mode,
                               parser, parser_cmd, config)
    texts: dict[Path, str] = {}
    for p, res in results.items():
        if res is None:
            continue
        _, migrated = res
        texts[p] = postprocess_migrated(migrated, p, config,
                                        module_rename_pairs)

    # Memoize the normalized text per path. A 4-member precision
    # family produces 3 pairs that all share the same canonical;
    # without this cache the canonical's _canonicalize_for_compare +
    # _strip_fortran_comments pipeline ran once per pair.
    normalized: dict[Path, str] = {}

    def _normalize(p: Path) -> str:
        n = normalized.get(p)
        if n is None:
            n = _canonicalize_for_compare(
                _strip_fortran_comments(texts[p], p.suffix))
            normalized[p] = n
        return n

    report: list[dict] = []
    for canonical, other in pairs:
        if canonical not in texts or other not in texts:
            continue
        n_can = _normalize(canonical)
        n_oth = _normalize(other)
        if n_can == n_oth:
            continue
        diff = _filter_precision_drift(
            n_oth.splitlines(), n_can.splitlines())
        if not diff:
            continue
        report.append({
            'target': target_filename(canonical.name, rename_map, target_mode),
            'canonical': canonical.name,
            'other': other.name,
            'diff': diff,
        })
    report.sort(key=lambda r: (r['other'], r['canonical']))
    if not apply_whitelist:
        return report
    return _filter_expected_divergences(report, config)
