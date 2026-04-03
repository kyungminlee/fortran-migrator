"""General-purpose C type migration engine.

Handles C source files that use the template-clone pattern (like BLACS):
type-specific files are near-identical clones differing only in C type
names and MPI datatype constants.

Migration is done by cloning files with mechanical text substitution.
No clang parser needed — C types are unambiguous single tokens.
"""

import re
import shutil
from pathlib import Path

from .prefix_classifier import target_prefix


# Default C type substitution rules.
# Each rule is (pattern, replacement_template).
# In replacement_template:
#   {REAL_TYPE}    → the target real type name (e.g., "QREAL", "EREAL")
#   {COMPLEX_TYPE} → the target complex type name (e.g., "QCOMPLEX")
#   {MPI_REAL}     → the target MPI real type (e.g., "MPI_QREAL")
#   {MPI_COMPLEX}  → the target MPI complex type (e.g., "MPI_QCOMPLEX")
#   {RP}           → real prefix char lowercase (e.g., "q", "e")
#   {CP}           → complex prefix char lowercase (e.g., "x", "y")
#   {RPU}          → real prefix char uppercase
#   {CPU}          → complex prefix char uppercase

REAL_CLONE_SUBS = [
    # Type names (run twice to catch adjacent matches)
    (r'(^|[^a-zA-Z_])double([^a-zA-Z_]|$)', r'\1{REAL_TYPE}\2'),
    (r'(^|[^a-zA-Z_])double([^a-zA-Z_]|$)', r'\1{REAL_TYPE}\2'),
    # MPI types
    (r'MPI_DOUBLE', 'MPI_{REAL_TYPE}'),
    # Function name prefixes
    (r'Cd([a-z])', r'C{RP}\1'),
    (r'BI_d([a-z])', r'BI_{RP}\1'),
]

COMPLEX_CLONE_SUBS = [
    # Complex struct types
    (r'(^|[^a-zA-Z_])DCOMPLEX([^a-zA-Z_]|$)', r'\1{COMPLEX_TYPE}\2'),
    (r'(^|[^a-zA-Z_])DCOMPLEX([^a-zA-Z_]|$)', r'\1{COMPLEX_TYPE}\2'),
    (r'(^|[^a-zA-Z_])SCOMPLEX([^a-zA-Z_]|$)', r'\1{COMPLEX_TYPE}\2'),
    (r'(^|[^a-zA-Z_])SCOMPLEX([^a-zA-Z_]|$)', r'\1{COMPLEX_TYPE}\2'),
    # Underlying real type
    (r'(^|[^a-zA-Z_])double([^a-zA-Z_]|$)', r'\1{REAL_TYPE}\2'),
    (r'(^|[^a-zA-Z_])double([^a-zA-Z_]|$)', r'\1{REAL_TYPE}\2'),
    # MPI types (order matters: DOUBLE_COMPLEX before DOUBLE)
    (r'MPI_DOUBLE_COMPLEX', 'MPI_{COMPLEX_TYPE}'),
    (r'MPI_DOUBLE', 'MPI_{REAL_TYPE}'),
    (r'MPI_COMPLEX([^a-zA-Z_0-9])', r'MPI_{COMPLEX_TYPE}\1'),
    # Function name prefixes
    (r'Cz([a-z])', r'C{CP}\1'),
    (r'BI_z([a-z])', r'BI_{CP}\1'),
    (r'BI_d([a-z])', r'BI_{RP}\1'),
]


def _build_sub_vars(kind: int) -> dict[str, str]:
    """Build template substitution variables for a given target KIND."""
    rp = target_prefix(kind, is_complex=False).lower()
    cp = target_prefix(kind, is_complex=True).lower()
    # Type names follow convention: Q/E for real prefix → QREAL/EREAL
    real_type = f'{rp.upper()}REAL'
    complex_type = f'{cp.upper()}COMPLEX'
    return {
        'REAL_TYPE': real_type,
        'COMPLEX_TYPE': complex_type,
        'MPI_REAL': f'MPI_{real_type}',
        'MPI_COMPLEX': f'MPI_{complex_type}',
        'RP': rp,
        'CP': cp,
        'RPU': rp.upper(),
        'CPU': cp.upper(),
    }


def clone_c_file(src_path: Path, dst_path: Path,
                 subs: list[tuple[str, str]],
                 template_vars: dict[str, str]) -> None:
    """Clone a C file with mechanical text substitutions."""
    text = src_path.read_text(errors='replace')

    for pattern, replacement in subs:
        # Expand template variables in replacement
        expanded = replacement
        for key, val in template_vars.items():
            expanded = expanded.replace(f'{{{key}}}', val)
        text = re.sub(pattern, expanded, text, flags=re.MULTILINE)

    dst_path.write_text(text)


def rename_c_file(name: str, old_prefix: str, new_prefix: str) -> str:
    """Rename a C file by replacing its prefix character."""
    stem = Path(name).stem
    ext = Path(name).suffix

    if stem.startswith(f'BI_{old_prefix}'):
        new_stem = f'BI_{new_prefix}' + stem[len(f'BI_{old_prefix}'):]
    elif stem.startswith(old_prefix):
        new_stem = new_prefix + stem[len(old_prefix):]
    else:
        return name

    return new_stem + ext


def migrate_c_directory(src_dir: Path, output_dir: Path,
                        kind: int, copy_originals: bool = True) -> dict:
    """Migrate a C source directory by cloning d→q/e and z→x/y variants.

    Returns a summary dict.
    """
    template_vars = _build_sub_vars(kind)
    rp = template_vars['RP']
    cp = template_vars['CP']

    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy all originals first
    if copy_originals:
        for f in sorted(src_dir.iterdir()):
            if f.suffix.lower() in ('.c', '.h'):
                shutil.copy2(f, output_dir / f.name)

    cloned = []

    # Clone d-variant → real-extended
    for f in sorted(src_dir.iterdir()):
        if f.suffix.lower() != '.c':
            continue
        stem = f.stem
        if stem.startswith('d') or stem.startswith('BI_d'):
            # Skip if it's not a type-specific file (check for 'd' prefix)
            if stem.startswith('BI_d') and stem[3] == 'd':
                pass  # BI_d* files
            elif stem.startswith('d') and not stem.startswith('BI_'):
                pass  # d* files
            else:
                continue

            new_name = rename_c_file(f.name, 'd', rp)
            if new_name == f.name:
                continue
            clone_c_file(f, output_dir / new_name,
                         REAL_CLONE_SUBS, template_vars)
            cloned.append(f'{f.name} → {new_name}')

    # Clone z-variant → complex-extended
    for f in sorted(src_dir.iterdir()):
        if f.suffix.lower() != '.c':
            continue
        stem = f.stem
        if stem.startswith('z') or stem.startswith('BI_z'):
            new_name = rename_c_file(f.name, 'z', cp)
            if new_name == f.name:
                continue
            clone_c_file(f, output_dir / new_name,
                         COMPLEX_CLONE_SUBS, template_vars)
            cloned.append(f'{f.name} → {new_name}')

    return {
        'cloned': cloned,
        'template_vars': template_vars,
    }
