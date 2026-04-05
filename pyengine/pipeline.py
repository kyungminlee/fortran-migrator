"""Migration pipeline — orchestrates the full migration for any library.

Usage:
    from pyengine.pipeline import run_migration
    run_migration(recipe_path, output_dir, target_kind=16)
"""

import os
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .config import RecipeConfig, load_recipe
from .symbol_scanner import scan_symbols
from .prefix_classifier import classify_symbols, build_rename_map
from .fortran_migrator import migrate_file, migrate_file_to_string, target_filename
from .c_migrator import migrate_c_directory

from tqdm import tqdm


def _strip_fortran_comments(text: str, ext: str) -> str:
    """Strip comments from Fortran source so convergence comparisons
    ignore commentary that differs between S/D or C/Z source halves.

    Fixed-form (.f/.for): a line is a full comment when column 1 is
    C, c, *, !, d, or D. Inline ! comments are stripped unless inside
    a string literal.

    Free-form (.f90/.F90/.f95): a line is a full comment when the
    first non-blank char is !. Inline ! comments are stripped the
    same way.

    Trailing whitespace and blank lines that result from stripping are
    removed so normalization is stable.
    """
    is_fixed = ext.lower() in ('.f', '.for')
    out_lines: list[str] = []
    for line in text.split('\n'):
        if is_fixed and line[:1] in ('C', 'c', '*', '!'):
            # Column-1 comment marker in fixed-form.
            continue
        if not is_fixed and line.lstrip().startswith('!'):
            continue
        # Strip inline ! comments (respect simple single/double quotes)
        stripped = _strip_inline_bang(line)
        if stripped.strip():
            out_lines.append(stripped.rstrip())

    # Join fixed-form continuation lines so sibling sources that wrap
    # at different columns still normalize to the same text. In
    # fixed-form, column 6 (non-blank, non-zero) indicates a
    # continuation of the previous line — merge it into its parent
    # after stripping the continuation marker and surrounding
    # whitespace.
    if is_fixed:
        joined: list[str] = []
        for ln in out_lines:
            if len(ln) > 5 and ln[5] not in (' ', '0', '\t') and joined:
                joined[-1] = joined[-1].rstrip() + ' ' + ln[6:].lstrip()
            else:
                joined.append(ln)
        out_lines = joined
    else:
        # Free-form: '&' at end of line continues onto next line.
        joined: list[str] = []
        pending = ''
        for ln in out_lines:
            if pending:
                ln = pending + ' ' + ln.lstrip()
                pending = ''
            stripped = ln.rstrip()
            if stripped.endswith('&'):
                pending = stripped[:-1].rstrip()
            else:
                joined.append(stripped)
        if pending:
            joined.append(pending)
        out_lines = joined
    return '\n'.join(out_lines)


def _strip_inline_bang(line: str) -> str:
    in_s = in_d = False
    for i, ch in enumerate(line):
        if ch == "'" and not in_d:
            in_s = not in_s
        elif ch == '"' and not in_s:
            in_d = not in_d
        elif ch == '!' and not in_s and not in_d:
            return line[:i]
    return line


def _canonicalize_for_compare(text: str) -> str:
    """Normalize text to ignore cosmetic precision-specific differences
    that remain after migration.

    Two sources of benign divergence are masked:

    * Numeric literal formatting — ``0.0``, ``0.0D0``, ``0.0E+0_16``,
      ``0.5_16`` etc. all represent the same value post-migration, but
      BLAS/LAPACK write them differently in S vs D sources. We
      canonicalize to a common form: replace ``[dD]`` exponent
      markers with ``E``, strip ``_N`` kind suffixes, and drop
      insignificant exponents (``E0``/``E+0``).

    * Prefix-dependent dummy-argument names — BLAS uses SA/SX/SY in
      S sources and DA/DX/DY in D sources (similarly CA/CX/CY vs
      ZA/ZX/ZY). These local names aren't in the rename map so they
      survive migration literally. We replace the leading S/D/C/Z
      of any identifier with a placeholder '@' so both halves look
      the same. (Applied symmetrically, so non-prefix identifiers
      like ``CALL``, ``COMPLEX``, ``SUBROUTINE`` are also mangled but
      mangled identically on both sides.)
    """
    # 1. d/D exponent → E in numeric literals
    text = re.sub(r'(\d)[dD]([+-]?\d)', r'\1E\2', text)
    # 2. Strip _N kind suffixes on literals (suffix after a digit)
    text = re.sub(r'(\d)_\d+\b', r'\1', text)
    # 3. Drop insignificant exponents (E0/E+0/e-0/E00 ...)
    text = re.sub(r'([0-9.])[eE][+-]?0+\b', r'\1', text)
    # 4. Canonicalize mantissa exponent marker to uppercase E
    text = re.sub(r'(\d)e([+-]?\d)', r'\1E\2', text)
    # 4b. Strip REAL(x, KIND=N) / CMPLX(x, KIND=N) wrappers around a
    #     bare identifier. Migrated S/C sources wrap every formerly
    #     single-precision REAL() cast with an explicit KIND=, while
    #     D/Z sources had no such cast in the original (they relied on
    #     implicit integer→double conversion). The two are semantically
    #     identical.
    text = re.sub(
        r'\b(?:REAL|CMPLX)\s*\(\s*(\w+)\s*,\s*KIND\s*=\s*\d+\s*\)',
        r'\1', text,
    )
    # 5. Canonicalize prefix-dependent identifiers: any run of leading
    #    S/D/C/Z characters on an identifier collapses to a single '@'.
    #    This handles both single-letter prefixes (SA vs DA) and
    #    double-prefix intrinsic names (CMPLX vs DCMPLX, CABS vs
    #    DCABS). Applied to both halves so the mapping is symmetric.
    text = re.sub(r'\b[sdczSDCZ]+(?=[A-Za-z])', '@', text)
    # 6. Collapse runs of whitespace to a single space — BLAS sources
    #    hand-align declaration lists with varying spacing between
    #    type keyword and argument list, and indent DO loop bodies
    #    differently across precision halves.
    text = re.sub(r'[ \t]+', ' ', text)
    # 7. Sort items in INTRINSIC / EXTERNAL lists. BLAS sources list
    #    intrinsics in arbitrary (often precision-specific) order —
    #    e.g., S source writes "INTRINSIC REAL, CONJG, MAX" while Z
    #    source writes "INTRINSIC CONJG, MAX, MIN, REAL". Sort both.
    def _sort_decl(match: re.Match) -> str:
        kw = match.group(1)
        items = [s.strip() for s in match.group(2).split(',')]
        return f'{kw} ' + ', '.join(sorted(items))
    text = re.sub(
        r'\b(INTRINSIC|EXTERNAL)\s+([A-Za-z_@][A-Za-z0-9_@,\s]*?)(?=$)',
        _sort_decl, text, flags=re.MULTILINE,
    )
    return text


def _strip_c_comments(text: str) -> str:
    """Strip // and /* */ comments plus blank/whitespace-only lines."""
    # Remove /* ... */ block comments (possibly multi-line)
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)
    # Remove // line comments
    text = re.sub(r'//[^\n]*', '', text)
    # Collapse blank/whitespace-only lines and trim trailing spaces
    lines = [ln.rstrip() for ln in text.split('\n')]
    lines = [ln for ln in lines if ln.strip()]
    return '\n'.join(lines)


def run_fortran_migration(config: RecipeConfig, rename_map: dict[str, str],
                          output_dir: Path, kind: int,
                          dry_run: bool = False,
                          classification=None) -> dict:
    """Run Fortran migration pipeline."""
    # Identify precision-independent symbols
    if classification is None:
        symbols = scan_symbols(config.source_dir, config.language,
                               config.extensions, config.library_path)
        classification = classify_symbols(symbols)
    independent = classification.independent

    migrated_count = 0
    copied_count = 0
    skipped: list[str] = []
    divergences: list[str] = []

    # Collect eligible source files
    src_files = sorted(
        p for p in config.source_dir.iterdir()
        if p.suffix.lower() in config.extensions
    )

    # Convergence buffer: first writer of each output name stores its
    # text; subsequent writers must agree or we record a divergence.
    # D/Z sources are preferred as the canonical text: when a pair
    # (SGEMM, DGEMM) both target QGEMM, DGEMM's migrated body is kept
    # and SGEMM's is only consulted for the equality check.
    def _is_double_source(stem: str) -> bool:
        return bool(stem) and stem[0].upper() in ('D', 'Z')

    # Process D/Z-sourced files first so their migration is the canonical
    # one on disk; S/C co-members are verified against them.
    src_files.sort(key=lambda p: (not _is_double_source(p.stem), p.name))

    canonical_text: dict[str, str] = {}
    canonical_normalized: dict[str, str] = {}
    canonical_source: dict[str, str] = {}

    # Partition: copy/skip/dry-run decisions stay in the main process;
    # only the Flang-bound migration is dispatched to workers.
    to_migrate: list[Path] = []
    for src_path in src_files:
        stem_upper = src_path.stem.upper()
        if stem_upper in config.skip_files:
            skipped.append(src_path.name)
            continue
        if dry_run:
            out_name = target_filename(src_path.name, rename_map)
            tqdm.write(f'  {src_path.name} → {out_name}')
            continue
        if stem_upper in config.copy_files or stem_upper in independent:
            (output_dir / src_path.name).write_text(
                src_path.read_text(errors='replace'))
            copied_count += 1
            continue
        to_migrate.append(src_path)

    if dry_run:
        return {
            'migrated': 0, 'copied': copied_count,
            'skipped': skipped, 'divergences': divergences,
        }

    # Parallel migration. Each worker runs flang + regex substitution,
    # which is the dominant cost for large libraries like LAPACK.
    # We reduce results in canonical-first order so D/Z output is
    # the one written to disk.
    workers = max(1, (os.cpu_count() or 4))
    results: dict[Path, tuple[str, str] | None] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(migrate_file_to_string, p, rename_map, kind): p
            for p in to_migrate
        }
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc='  Migrating', unit='file',
                        mininterval=1.0, miniters=10):
            p = futures[fut]
            results[p] = fut.result()

    # Reduce in deterministic D/Z-first order.
    for src_path in to_migrate:
        result = results.get(src_path)
        if result is None:
            skipped.append(src_path.name)
            continue
        out_name, migrated = result
        normalized = _canonicalize_for_compare(
            _strip_fortran_comments(migrated, src_path.suffix)
        )

        prior = canonical_normalized.get(out_name)
        if prior is None:
            # First (canonical) writer for this target name
            (output_dir / out_name).write_text(migrated)
            canonical_text[out_name] = migrated
            canonical_normalized[out_name] = normalized
            canonical_source[out_name] = src_path.name
            migrated_count += 1
        elif prior == normalized:
            # Convergence: co-family member produced identical code
            # (comments may differ and are ignored).
            pass
        else:
            # Divergence: co-family members disagree. Keep the canonical
            # (D/Z) version on disk and record the mismatch.
            divergences.append(
                f'{src_path.name} vs {canonical_source[out_name]} '
                f'→ {out_name}'
            )

    return {
        'migrated': migrated_count,
        'copied': copied_count,
        'skipped': skipped,
        'divergences': divergences,
    }


def run_c_migration(config: RecipeConfig, output_dir: Path,
                    kind: int, dry_run: bool = False,
                    classification=None,
                    rename_map: dict[str, str] | None = None) -> dict:
    """Run C migration pipeline (clone-and-substitute).

    When `classification` and `rename_map` are supplied, the generic
    rename-map-driven migration is used (for ScaLAPACK-style libraries
    like PBLAS). Otherwise the BLACS-specific path is used.
    """
    if dry_run:
        print('  (dry-run for C migration not yet implemented)')
        return {'cloned': [], 'template_vars': {}}

    result = migrate_c_directory(
        config.source_dir, output_dir, kind,
        copy_originals=config.copy_all_originals,
        classification=classification,
        rename_map=rename_map,
    )
    return result


def run_migration(recipe_path: Path, output_dir: Path,
                  target_kind: int = 16, dry_run: bool = False,
                  project_root: Path | None = None) -> dict:
    """Run the full migration pipeline for a library.

    Args:
        recipe_path: Path to the YAML recipe file.
        output_dir: Where to write migrated files.
        target_kind: 10 or 16.
        dry_run: If True, show what would be done without writing.
        project_root: Project root for resolving relative paths.

    Returns:
        Summary dict with migration statistics.
    """
    config = load_recipe(recipe_path, project_root)
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f'Library:     {config.library}')
    print(f'Language:    {config.language}')
    print(f'Source:      {config.source_dir}')
    print(f'Target KIND: {target_kind}')
    print(f'Prefix:      {config.prefix_style}')
    print()

    # Scan symbols and classify precision families
    print('Scanning symbols...')
    symbols = scan_symbols(
        config.source_dir, config.language,
        config.extensions, config.library_path
    )
    own_count = len(symbols)

    # Merge symbols from dependency libraries
    for dep_path in config.depends:
        dep_config = load_recipe(dep_path, project_root)
        dep_symbols = scan_symbols(
            dep_config.source_dir, dep_config.language,
            dep_config.extensions, dep_config.library_path
        )
        symbols |= dep_symbols
        # Recursively load transitive dependencies
        for transitive in dep_config.depends:
            trans_config = load_recipe(transitive, project_root)
            trans_symbols = scan_symbols(
                trans_config.source_dir, trans_config.language,
                trans_config.extensions, trans_config.library_path
            )
            symbols |= trans_symbols

    # Scan extra symbol directories (for external dependencies like PBLAS/TOOLS)
    # These are scanned recursively to handle subdirectory structures.
    # Both Fortran and C files are scanned for symbols.
    for extra_dir in config.extra_symbol_dirs:
        if extra_dir.is_dir():
            for lang, exts in [('fortran', ['.f', '.f90', '.F90']),
                               ('c', ['.c'])]:
                extra_syms = scan_symbols(extra_dir, lang, exts)
                symbols |= extra_syms
            # Also scan subdirectories
            for sub in sorted(extra_dir.iterdir()):
                if sub.is_dir():
                    for lang, exts in [('fortran', ['.f', '.f90', '.F90']),
                                       ('c', ['.c'])]:
                        sub_syms = scan_symbols(sub, lang, exts)
                        symbols |= sub_syms

    classification = classify_symbols(symbols)
    rename_map = classification.build_rename_map(target_kind)
    print(f'  {own_count} own symbols + {len(symbols) - own_count} from dependencies')
    print(f'  {len(classification.families)} precision families')
    print(f'  {len(classification.independent)} independent symbols')
    print(f'  {len(rename_map)} renames computed')

    # Dispatch to language-specific migrator
    print(f'\nMigrating to KIND={target_kind}...')

    if config.language == 'fortran':
        result = run_fortran_migration(
            config, rename_map, output_dir, target_kind, dry_run,
            classification=classification
        )
        if not dry_run:
            print(f'\n  Migrated: {result["migrated"]} files')
            print(f'  Copied:   {result["copied"]} files (precision-independent)')
            if result['skipped']:
                print(f'  Skipped:  {len(result["skipped"])} files')
            if result.get('divergences'):
                print(f'  Divergences: {len(result["divergences"])} '
                      f'(co-family members produced differing output)')
                for d in result['divergences'][:10]:
                    print(f'    {d}')
                if len(result['divergences']) > 10:
                    print(f'    ... ({len(result["divergences"]) - 10} more)')
    elif config.language == 'c':
        # ScaLAPACK-style C libraries (PBLAS) use rename-map-driven
        # cloning; BLACS-style (prefix 'direct') keeps the legacy path.
        if config.prefix_style == 'scalapack':
            result = run_c_migration(
                config, output_dir, target_kind, dry_run,
                classification=classification, rename_map=rename_map,
            )
        else:
            result = run_c_migration(config, output_dir, target_kind, dry_run)
        if not dry_run:
            print(f'\n  Cloned: {len(result["cloned"])} files')
            if result.get('divergences'):
                print(f'  Divergences: {len(result["divergences"])} '
                      f'(co-family members produced differing output)')
                for d in result['divergences'][:10]:
                    print(f'    {d}')
                if len(result['divergences']) > 10:
                    print(f'    ... ({len(result["divergences"]) - 10} more)')
    else:
        raise ValueError(f'Unsupported language: {config.language}')

    return result
