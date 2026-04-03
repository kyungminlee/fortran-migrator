"""CLI entry point: python -m pyengine <recipe.yaml> <output-dir> [--kind 16]"""

import argparse
import sys
from pathlib import Path

from .pipeline import run_migration


def main():
    parser = argparse.ArgumentParser(
        prog='pyengine',
        description='General-purpose type migration for numerical libraries'
    )
    parser.add_argument('recipe', type=Path,
                        help='Path to recipe YAML file')
    parser.add_argument('output_dir', type=Path,
                        help='Output directory for migrated files')
    parser.add_argument('--kind', type=int, default=16, choices=[10, 16],
                        help='Target floating-point KIND (default: 16)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show planned changes without writing')
    parser.add_argument('--project-root', type=Path, default=None,
                        help='Project root for resolving relative paths '
                             '(default: recipe parent\'s parent)')
    args = parser.parse_args()

    if not args.recipe.exists():
        print(f'Error: recipe not found: {args.recipe}', file=sys.stderr)
        sys.exit(1)

    run_migration(
        recipe_path=args.recipe,
        output_dir=args.output_dir,
        target_kind=args.kind,
        dry_run=args.dry_run,
        project_root=args.project_root,
    )


if __name__ == '__main__':
    main()
