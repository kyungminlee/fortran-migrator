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


def _strip_real_cmplx_casts(text: str) -> str:
    """Strip ``REAL(expr, KIND=N)`` and ``CMPLX(expr, KIND=N)`` wrappers,
    replacing them with just ``expr``. Single-argument ``REAL(expr)``
    (integer/implicit-real promotion) is also stripped, but
    ``CMPLX(expr)`` is NOT — it materially changes the result type
    from real to complex and the two halves genuinely may differ in
    how they construct complex literals."""
    pattern = re.compile(r'\b(?:REAL|CMPLX)\s*\(', re.IGNORECASE)
    out = []
    i = 0
    while i < len(text):
        m = pattern.search(text, i)
        if not m:
            out.append(text[i:])
            break
        out.append(text[i:m.start()])
        # Find matching ')'
        depth = 1
        j = m.end()
        while j < len(text) and depth > 0:
            c = text[j]
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth != 0:
            # Unmatched — give up on this match
            out.append(text[m.start():m.end()])
            i = m.end()
            continue
        inner = text[m.end():j]
        # Recurse into the inner expression so nested casts are
        # stripped even when the outer call doesn't match the KIND
        # shape we rewrite.
        inner = _strip_real_cmplx_casts(inner)
        parts = _split_top_level_comma(inner)
        is_real_call = text[m.start():m.end()].upper().lstrip().startswith('REAL')
        # Skip declaration-form ``REAL(KIND=N)`` — that's a type
        # specifier, not a cast.
        is_kind_decl = (len(parts) == 1 and
                        re.fullmatch(r'\s*KIND\s*=\s*\d+\s*', parts[0]))
        # Strip trailing ``KIND=N`` from any multi-arg call. This
        # covers ``CMPLX(re, im, KIND=16)`` as well as
        # ``REAL(x, KIND=16)``.
        is_kinded_multi = (len(parts) >= 2 and
                           re.fullmatch(r'\s*KIND\s*=\s*\d+\s*', parts[-1]))
        if is_kinded_multi:
            # Drop the KIND= tail argument but keep the call.
            call_head = text[m.start():m.end()]
            new_inner = ','.join(p.strip() for p in parts[:-1])
            if len(parts) == 2:
                # 1-data-arg form: reduce to bare identifier/expr.
                inner_expr = parts[0].strip()
                if _has_top_level_operator(inner_expr):
                    out.append('(' + inner_expr + ')')
                else:
                    out.append(inner_expr)
            else:
                # Multi-arg CMPLX(a, b, KIND=N) → CMPLX(a, b)
                out.append(call_head + new_inner + ')')
            i = j + 1
            continue
        if not is_kind_decl and (
                (len(parts) == 1 and is_real_call)):
            # Wrap the stripped expression in parens ONLY when the
            # inner contains a top-level additive operator — the
            # surrounding context's precedence then matches the
            # unwrapped side (``REAL(N-1, KIND=16)*T`` → ``(N-1)*T``
            # rather than ``N-1*T``). Single atoms keep their natural
            # form so ``REAL(X, KIND=16)`` → ``X``.
            #
            # Single-arg REAL(X)/CMPLX(X) is also stripped: it is an
            # integer→real promotion that the Fortran compiler would
            # insert implicitly when X is used in a real context, so
            # ``REAL(MAXITR)*UNFL`` is semantically identical to
            # ``MAXITR*UNFL`` after migration.
            inner_expr = parts[0].strip()
            if _has_top_level_operator(inner_expr):
                out.append('(' + inner_expr + ')')
            else:
                out.append(inner_expr)
        else:
            # Not the shape we want — keep original call, but with
            # inner expression already stripped recursively.
            out.append(text[m.start():m.end()] + inner + ')')
        i = j + 1
    return ''.join(out)


def _has_top_level_operator(s: str) -> bool:
    """True if ``s`` contains a top-level ``+`` or ``-`` outside nested
    parens/brackets (ignoring unary signs and exponent markers).

    We only consider additive operators here, not ``*`` / ``/``:
    when the stripped expression is used as a factor in an outer
    multiplication (the usual case for REAL(int_expr, KIND=16)*x),
    a pure product like ``MAXITR*Q*Q`` is associative with the
    outer ``*`` so parens are unnecessary. Only ``a-b`` or ``a+b``
    would change precedence under multiplication.
    """
    depth = 0
    for i, c in enumerate(s):
        if c in '([':
            depth += 1
        elif c in ')]':
            depth -= 1
        elif depth == 0 and c in '+-' and i > 0:
            # Skip if preceded by 'E' / 'D' (exponent sign) or an
            # operator (meaning it's a unary sign).
            prev = s[i - 1]
            if prev in 'eEdD' and i >= 2 and s[i - 2].isdigit():
                continue
            if prev in '+-*/=(,':
                continue
            return True
    return False


def _split_top_level_comma(s: str) -> list[str]:
    """Split on commas that are not inside nested parens or brackets."""
    out: list[str] = []
    depth = 0
    start = 0
    for i, c in enumerate(s):
        if c in '([':
            depth += 1
        elif c in ')]':
            depth -= 1
        elif c == ',' and depth == 0:
            out.append(s[start:i])
            start = i + 1
    out.append(s[start:])
    return out


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
    # 4b. Strip REAL(expr, KIND=N) / CMPLX(expr, KIND=N) wrappers.
    #     Migrated S/C sources wrap every formerly single-precision
    #     REAL() cast with an explicit KIND=, while D/Z sources had no
    #     such cast in the original (they relied on implicit
    #     integer→real conversion). The two are semantically
    #     identical. We match balanced parens because the inner
    #     expression may contain nested calls like ``REAL(F(I)*G, KIND=16)``.
    text = _strip_real_cmplx_casts(text)
    # 4c. Canonicalize numeric literal forms. ``1.0`` / ``1.`` / ``1``
    #     all denote the same value; LAPACK writes them inconsistently
    #     between S/C and D/Z halves (``1.0,0.0`` complex zero vs
    #     ``1.,0.``). Drop fractional parts whose numeric value is
    #     zero so all three forms collapse to ``1``.
    text = re.sub(r'(\d)\.0*(?![\deEdD_])', r'\1', text)
    # 4d. Strip ``(KIND=N)`` suffix on REAL/COMPLEX type-spec in
    #     F90 declarations. Some LAPACK F90 files use the bare
    #     ``COMPLEX, INTENT(...)`` form while the migrator adds
    #     ``(KIND=16)`` to others — normalize to bare form. Must run
    #     BEFORE the @-collapse below, which would otherwise turn
    #     ``COMPLEX`` into ``@OMPLEX`` and prevent this regex from
    #     matching.
    text = re.sub(
        r'\b(REAL|COMPLEX)\s*\(\s*KIND\s*=\s*\d+\s*\)',
        r'\1', text, flags=re.IGNORECASE)
    # 5. Canonicalize prefix-dependent identifiers: any run of leading
    #    S/D/C/Z characters on an identifier collapses to a single '@'.
    #    This handles both single-letter prefixes (SA vs DA) and
    #    double-prefix intrinsic names (CMPLX vs DCMPLX, CABS vs
    #    DCABS). Applied to both halves so the mapping is symmetric.
    text = re.sub(
        r'(?<![A-Za-z0-9])[sdczSDCZ]+(?=[A-Za-z])', '@', text,
    )
    # 5a. iso_fortran_env kind imports: ``only: real32`` vs
    #     ``only: real64`` survive as an unused import (after the
    #     migrator rewrites ``WP = real32`` → ``WP = 16``). Normalize
    #     the kind-suffix digits so both halves compare equal.
    text = re.sub(r'\breal(32|64|128)\b', 'realK', text, flags=re.IGNORECASE)
    # 5b. Uppercase everything. Fortran is case-insensitive, and S vs
    #     D sources are not consistent about casing keywords
    #     (``then`` vs ``THEN``, ``IF`` vs ``if``). Uppercasing both
    #     halves makes them compare equal.
    text = text.upper()
    # 5c. Merge ``END IF`` → ``ENDIF`` etc. LAPACK writes these block
    #     terminators with or without a space arbitrarily between
    #     precision halves.
    text = re.sub(r'\bEND\s+(IF|DO|SELECT|WHERE|FORALL|SUBROUTINE|FUNCTION|MODULE|INTERFACE|PROGRAM|TYPE|BLOCK|ASSOCIATE)\b',
                  r'END\1', text)
    # 5d. Strip bare ``IMPLICIT NONE`` (or any other IMPLICIT spec).
    #     Some S/D halves have it and others don't; it has no bearing
    #     on the migrated numerics.
    text = re.sub(r'^\s*IMPLICIT\s+\S.*$', '', text, flags=re.MULTILINE)
    # 6. Collapse runs of whitespace to a single space — BLAS sources
    #    hand-align declaration lists with varying spacing between
    #    type keyword and argument list, and indent DO loop bodies
    #    differently across precision halves.
    text = re.sub(r'[ \t]+', ' ', text)
    # 6b. Strip whitespace adjacent to arithmetic / relational
    #     operators. LAPACK's S and D sources format the same
    #     expression with different spacing (``NZ*EPS`` vs
    #     ``NZ* EPS``, ``I + 1`` vs ``I+1``). We squeeze both sides
    #     so they compare equal.
    text = re.sub(r'\s*([*+\-/=])\s*', r'\1', text)
    # 6c. Strip whitespace adjacent to parens. LAPACK formats
    #     argument lists with varying spacing (``F( X )`` vs
    #     ``F(X)``, ``IF (X)`` vs ``IF(X)``).
    text = re.sub(r'\s*\(\s*', '(', text)
    text = re.sub(r'\s*\)', ')', text)
    # 6d. Strip whitespace after commas. LAPACK formats complex
    #     literals and argument lists as ``(1.0, 0.0)`` vs
    #     ``(1.0,0.0)`` inconsistently between S/C and D/Z halves.
    text = re.sub(r',\s+', ',', text)
    # 6e. Collapse doubled parens ``((X))`` → ``(X)``. LAPACK
    #     occasionally wraps a boolean / logical expression in an
    #     extra redundant paren on one half.
    prev = None
    while prev != text:
        prev = text
        text = re.sub(r'\(\(([^()]*)\)\)', r'(\1)', text)
    # 7. Sort items in declaration lists whose order isn't
    #    semantically meaningful. BLAS/LAPACK sources list
    #    intrinsics, externals, and type-declared variables in
    #    arbitrary (often precision-specific) order — e.g., S source
    #    writes "INTRINSIC REAL, CONJG, MAX" while Z source writes
    #    "INTRINSIC CONJG, MAX, MIN, REAL".
    def _sort_decl(match: re.Match) -> str:
        kw = match.group(1)
        items = [s.strip() for s in match.group(2).split(',')]
        return f'{kw} ' + ', '.join(sorted(items))
    text = re.sub(
        r'\b(INTRINSIC|EXTERNAL|INTEGER|LOGICAL|COMPLEX|REAL)\s+'
        r'([A-Za-z_@][A-Za-z0-9_@,\s]*?)(?=$)',
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


def run_divergence_report(recipe_path: Path, target_kind: int = 16,
                           project_root: Path | None = None) -> list[dict]:
    """Migrate every co-family source pair in-memory and return the
    normalized diff for each pair whose members disagree.

    Returns a list of ``{'target', 'canonical', 'other', 'diff'}``
    dicts, sorted by target filename. ``diff`` is the list of +/-
    lines (without context) from the unified diff of the two
    canonicalized texts.
    """
    import difflib
    config = load_recipe(recipe_path, project_root)

    # Scan this recipe and its dependencies.
    symbols = scan_symbols(
        config.source_dir, config.language,
        config.extensions, config.library_path,
    )
    for dep_path in config.depends:
        dep_cfg = load_recipe(dep_path, project_root)
        symbols |= scan_symbols(
            dep_cfg.source_dir, dep_cfg.language,
            dep_cfg.extensions, dep_cfg.library_path,
        )
    for extra_dir in config.extra_symbol_dirs:
        if extra_dir.is_dir():
            for lang, exts in [('fortran', ['.f', '.f90', '.F90']),
                               ('c', ['.c'])]:
                symbols |= scan_symbols(extra_dir, lang, exts)
            for sub in sorted(extra_dir.iterdir()):
                if sub.is_dir():
                    for lang, exts in [('fortran', ['.f', '.f90', '.F90']),
                                       ('c', ['.c'])]:
                        symbols |= scan_symbols(sub, lang, exts)

    classification = classify_symbols(symbols)
    rename_map = classification.build_rename_map(target_kind)

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
        by_target.setdefault(target_filename(p.name, rename_map), []).append(p)

    # Migrate every member of every multi-member group in parallel.
    pairs: list[tuple[Path, Path]] = []
    for tn, members in by_target.items():
        if len(members) < 2:
            continue
        members.sort(key=lambda p: (p.stem[0].upper() not in ('D', 'Z'), p.name))
        canonical = members[0]
        for other in members[1:]:
            pairs.append((canonical, other))

    all_paths = {p for pair in pairs for p in pair}
    texts: dict[Path, str] = {}
    workers = max(1, (os.cpu_count() or 4))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(migrate_file_to_string, p, rename_map, target_kind): p
            for p in all_paths
        }
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc='  Migrating', unit='file',
                        mininterval=1.0, miniters=10):
            p = futures[fut]
            try:
                _, migrated = fut.result()
                texts[p] = migrated
            except Exception:
                pass

    report: list[dict] = []
    for canonical, other in pairs:
        if canonical not in texts or other not in texts:
            continue
        n_can = _canonicalize_for_compare(
            _strip_fortran_comments(texts[canonical], canonical.suffix))
        n_oth = _canonicalize_for_compare(
            _strip_fortran_comments(texts[other], other.suffix))
        if n_can == n_oth:
            continue
        diff = [
            line for line in difflib.unified_diff(
                n_oth.splitlines(), n_can.splitlines(),
                lineterm='', n=0,
            )
            if line.startswith(('-', '+'))
            and not line.startswith(('---', '+++'))
        ]
        report.append({
            'target': target_filename(canonical.name, rename_map),
            'canonical': canonical.name,
            'other': other.name,
            'diff': diff,
        })
    report.sort(key=lambda r: (r['other'], r['canonical']))
    return report


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
