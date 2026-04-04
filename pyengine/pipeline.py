"""Migration pipeline — orchestrates the full migration for any library.

Usage:
    from pyengine.pipeline import run_migration
    run_migration(recipe_path, output_dir, target_kind=16)
"""

from pathlib import Path

from .config import RecipeConfig, load_recipe
from .symbol_scanner import scan_symbols
from .prefix_classifier import classify_symbols, build_rename_map
from .fortran_migrator import migrate_file, target_filename
from .c_migrator import migrate_c_directory

from tqdm import tqdm


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

    # Collect eligible source files
    src_files = sorted(
        p for p in config.source_dir.iterdir()
        if p.suffix.lower() in config.extensions
    )

    for src_path in tqdm(src_files, desc='  Migrating', unit='file',
                          mininterval=1.0, miniters=10):
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

        out_name = migrate_file(src_path, output_dir, rename_map, kind)
        if out_name:
            migrated_count += 1
        else:
            skipped.append(src_path.name)

    return {
        'migrated': migrated_count,
        'copied': copied_count,
        'skipped': skipped,
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
    else:
        raise ValueError(f'Unsupported language: {config.language}')

    return result
