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
    (r'MPI_DOUBLE', '{MPI_REAL}'),
    # Function name prefixes (allow uppercase after prefix for BI_dMPI_* etc.)
    (r'Cd([a-z])', r'C{RP}\1'),
    (r'BI_d([a-zA-Z])', r'BI_{RP}\1'),
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
    (r'MPI_DOUBLE_COMPLEX', '{MPI_COMPLEX}'),
    (r'MPI_DOUBLE', '{MPI_REAL}'),
    (r'MPI_COMPLEX([^a-zA-Z_0-9])', r'{MPI_COMPLEX}\1'),
    # Function name prefixes (allow uppercase after prefix for BI_zMPI_* etc.)
    (r'Cz([a-z])', r'C{CP}\1'),
    (r'BI_z([a-zA-Z])', r'BI_{CP}\1'),
    (r'BI_d([a-zA-Z])', r'BI_{RP}\1'),
]


# Standard MPI datatype names for each target KIND.
# KIND=16 (quad precision): MPI_REAL16 / MPI_COMPLEX32 — requires MPI
#   implementation with quad-precision support.
# KIND=10 (extended precision): maps to long double on x86.
_MPI_REAL_TYPE = {
    16: 'MPI_REAL16',
    10: 'MPI_LONG_DOUBLE',
}
_MPI_COMPLEX_TYPE = {
    16: 'MPI_COMPLEX32',
    10: 'MPI_C_LONG_DOUBLE_COMPLEX',
}


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
        'MPI_REAL': _MPI_REAL_TYPE[kind],
        'MPI_COMPLEX': _MPI_COMPLEX_TYPE[kind],
        'RP': rp,
        'CP': cp,
        'RPU': rp.upper(),
        'CPU': cp.upper(),
    }


def clone_c_file(src_path: Path, dst_path: Path,
                 subs: list[tuple[str, str]],
                 template_vars: dict[str, str],
                 routine_renames: list[tuple[str, str]] | None = None,
                 ) -> None:
    """Clone a C file with mechanical text substitutions.

    routine_renames is a list of (old_name, new_name) pairs for literal
    routine name replacements (both lowercase and uppercase are applied).
    """
    text = src_path.read_text(errors='replace')

    for pattern, replacement in subs:
        # Expand template variables in replacement
        expanded = replacement
        for key, val in template_vars.items():
            expanded = expanded.replace(f'{{{key}}}', val)
        text = re.sub(pattern, expanded, text, flags=re.MULTILINE)

    # Apply routine name renames (lowercase and uppercase)
    if routine_renames:
        for old_name, new_name in routine_renames:
            text = text.replace(old_name, new_name)
            text = text.replace(old_name.upper(), new_name.upper())

    dst_path.write_text(text)


def _routine_renames(old_stem: str, new_stem: str) -> list[tuple[str, str]]:
    """Derive routine name renames from source/target file stems.

    For user-facing files like 'dgesd2d_' → 'qgesd2d_', strips the
    trailing underscore to get the routine base name and returns rename
    pairs.  For BI_-prefixed files the function names are already handled
    by the regex rules, so this returns an empty list.
    """
    if old_stem.startswith('BI_'):
        return []
    # Strip trailing underscore if present (Fortran naming convention)
    old_routine = old_stem.rstrip('_')
    new_routine = new_stem.rstrip('_')
    if old_routine == new_routine:
        return []
    return [(old_routine, new_routine)]


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
            renames = _routine_renames(stem, Path(new_name).stem)
            clone_c_file(f, output_dir / new_name,
                         REAL_CLONE_SUBS, template_vars, renames)
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
            renames = _routine_renames(stem, Path(new_name).stem)
            clone_c_file(f, output_dir / new_name,
                         COMPLEX_CLONE_SUBS, template_vars, renames)
            cloned.append(f'{f.name} → {new_name}')

    # Generate MPI requirement check for KIND=16
    if kind == 16:
        _generate_mpi_real16_check(output_dir)

    return {
        'cloned': cloned,
        'template_vars': template_vars,
    }


def _generate_mpi_real16_check(output_dir: Path) -> None:
    """Generate a CMake module that verifies MPI_REAL16 support."""
    cmake = """\
# CheckMpiReal16.cmake — Verify that the MPI implementation provides
# MPI_REAL16 and MPI_COMPLEX32, which are required by the migrated
# quad-precision BLACS routines.
#
# Usage:
#   include(CheckMpiReal16)
#   check_mpi_real16()          # FATAL_ERROR if unsupported
#
# After a successful check the cache variable HAVE_MPI_REAL16 is set.

include(CheckCSourceCompiles)

function(check_mpi_real16)
    find_package(MPI REQUIRED COMPONENTS C)

    set(CMAKE_REQUIRED_INCLUDES  ${MPI_C_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS     ${MPI_C_COMPILE_FLAGS})

    check_c_source_compiles("
        #include <mpi.h>
        int main(void) {
            MPI_Datatype dt_real = MPI_REAL16;
            MPI_Datatype dt_cplx = MPI_COMPLEX32;
            (void)dt_real; (void)dt_cplx;
            return 0;
        }
    " HAVE_MPI_REAL16)

    if(NOT HAVE_MPI_REAL16)
        message(FATAL_ERROR
            "The MPI implementation does not provide MPI_REAL16 / MPI_COMPLEX32. "
            "Quad-precision BLACS requires an MPI library built with 128-bit real "
            "support (e.g. OpenMPI or MPICH compiled with __float128 enabled).")
    endif()

    set(HAVE_MPI_REAL16 ${HAVE_MPI_REAL16} PARENT_SCOPE)
endfunction()
"""
    (output_dir / 'CheckMpiReal16.cmake').write_text(cmake)
