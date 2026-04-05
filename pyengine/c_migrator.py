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

from .prefix_classifier import target_prefix, SymbolClassification


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
                        kind: int, copy_originals: bool = True,
                        classification: SymbolClassification | None = None,
                        rename_map: dict[str, str] | None = None) -> dict:
    """Migrate a C source directory by cloning real/complex variants.

    Two modes:

    - **Generic** (when `classification` and `rename_map` are supplied):
      clones each D/Z member of every precision family found in
      `classification`, applying `rename_map` for all in-file routine
      name substitutions plus C-type upgrades (``double`` →
      ``{REAL_TYPE}``). Used for ScaLAPACK-style libraries like PBLAS
      whose entry points are Fortran-callable names such as
      ``pdgemm_``, ``pzhemm_``, ``pdznrm2_``.

    - **BLACS** (when neither is supplied): hardcoded ``d → q`` / ``z →
      x`` file renames with BLACS-specific ``Cd*``/``BI_d*`` routine
      patterns, MPI type substitutions, ``Bdef.h`` patching, and an
      MPI_REAL16 check module for KIND=16. Preserved for backward
      compatibility.

    Returns a summary dict.
    """
    if classification is not None and rename_map is not None:
        return _migrate_generic_c_directory(
            src_dir, output_dir, kind, copy_originals,
            classification, rename_map,
        )
    return _migrate_blacs_c_directory(src_dir, output_dir, kind, copy_originals)


def _build_rename_regex(rename_map: dict[str, str]) -> tuple[re.Pattern, dict[str, str]]:
    """Build a single-pass regex that renames routine names in C text.

    For each (OLD, NEW) in rename_map we emit four lookup entries:
    uppercase bare, lowercase bare, uppercase with trailing underscore,
    lowercase with trailing underscore. Word boundaries prevent false
    matches inside longer identifiers (e.g. ``DGER`` inside ``PDGER``).
    """
    combined: dict[str, str] = {}
    for old, new in rename_map.items():
        combined[old] = new
        combined[old.lower()] = new.lower()
        combined[old + '_'] = new + '_'
        combined[old.lower() + '_'] = new.lower() + '_'
    # Longest keys first so the alternation prefers the underscore forms
    keys_sorted = sorted(combined.keys(), key=len, reverse=True)
    pattern = re.compile(
        r'\b(' + '|'.join(re.escape(k) for k in keys_sorted) + r')\b'
    )
    return pattern, combined


def _apply_c_type_subs(text: str, template_vars: dict[str, str]) -> str:
    """Upgrade C type names used by precision-specific source files."""
    text = re.sub(r'\bdouble\b', template_vars['REAL_TYPE'], text)
    text = re.sub(r'\bfloat\b', template_vars['REAL_TYPE'], text)
    text = re.sub(r'\bDCOMPLEX\b', template_vars['COMPLEX_TYPE'], text)
    text = re.sub(r'\bSCOMPLEX\b', template_vars['COMPLEX_TYPE'], text)
    return text


def _migrate_generic_c_directory(src_dir: Path, output_dir: Path,
                                 kind: int, copy_originals: bool,
                                 classification: SymbolClassification,
                                 rename_map: dict[str, str]) -> dict:
    """Rename-map-driven C migration for ScaLAPACK-style libraries.

    A file ``foo_.c`` is cloned iff its routine ``FOO`` is a D- or
    Z-precision member of some family in ``classification``. The clone
    is written to ``<target>.c`` where ``<target>`` is the lowercase of
    ``rename_map[FOO]``, with all in-file routine names rewritten via
    ``rename_map`` and C types upgraded to the target precision.
    """
    template_vars = _build_sub_vars(kind)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy all originals first (keeps S/D/C/Z entry points available)
    if copy_originals:
        for f in sorted(src_dir.iterdir()):
            if f.suffix.lower() in ('.c', '.h'):
                shutil.copy2(f, output_dir / f.name)

    rename_pattern, combined_map = _build_rename_regex(rename_map)

    def _rename(text: str) -> str:
        return rename_pattern.sub(lambda m: combined_map[m.group(0)], text)

    # Process D/Z-sourced files first so they become the canonical
    # output; S/C co-family members are verified against them.
    def _is_double_key(routine_upper: str) -> bool:
        fam = classification.get_family(routine_upper)
        if fam is None:
            return False
        key = next(
            (k for k, v in fam.members.items() if v == routine_upper),
            '',
        )
        return bool(key) and key[0] in ('D', 'Z')

    entries = sorted(
        (f for f in src_dir.iterdir() if f.suffix.lower() == '.c'),
        key=lambda f: (
            not _is_double_key(
                (f.stem[:-1] if f.stem.endswith('_') else f.stem).upper()
            ),
            f.name,
        ),
    )

    def _normalize_for_compare(s: str) -> str:
        # Strip block + line comments
        s = re.sub(r'/\*.*?\*/', '', s, flags=re.DOTALL)
        s = re.sub(r'//[^\n]*', '', s)
        # Canonicalize prefix-dependent identifiers: leading s/d/c/z
        # becomes '@' so sibling C sources that differ only in the
        # precision prefix of local names collapse together.
        s = re.sub(r'\b[sdczSDCZ]+(?=[A-Za-z])', '@', s)
        # Collapse all whitespace so column-aligned declarations like
        # ``float          * ALPHA;`` and ``double         * ALPHA;``
        # compare equal after the type substitution.
        lines = [re.sub(r'\s+', ' ', ln).strip() for ln in s.split('\n')]
        return '\n'.join(ln for ln in lines if ln)

    cloned: list[str] = []
    divergences: list[str] = []
    canonical_normalized: dict[str, str] = {}
    canonical_source: dict[str, str] = {}

    for f in entries:
        has_underscore = f.stem.endswith('_')
        routine = f.stem[:-1] if has_underscore else f.stem
        upper_routine = routine.upper()

        if upper_routine not in rename_map:
            continue
        if classification.get_family(upper_routine) is None:
            continue

        target_upper = rename_map[upper_routine]
        target_lower = target_upper.lower()
        new_stem = target_lower + ('_' if has_underscore else '')
        new_name = new_stem + f.suffix
        new_path = output_dir / new_name

        text = f.read_text(errors='replace')
        # Apply renames first, then type upgrades — the two domains
        # don't overlap but this order preserves identifier names that
        # happen to coincide with generic type keywords.
        text = _rename(text)
        text = _apply_c_type_subs(text, template_vars)

        normalized = _normalize_for_compare(text)
        prior = canonical_normalized.get(new_name)
        if prior is None:
            new_path.write_text(text)
            canonical_normalized[new_name] = normalized
            canonical_source[new_name] = f.name
            cloned.append(f'{f.name} → {new_name}')
        elif prior == normalized:
            pass  # convergence (ignoring comment differences)
        else:
            divergences.append(
                f'{f.name} vs {canonical_source[new_name]} → {new_name}'
            )

    return {
        'cloned': cloned,
        'divergences': divergences,
        'template_vars': template_vars,
    }


def _migrate_blacs_c_directory(src_dir: Path, output_dir: Path,
                               kind: int, copy_originals: bool) -> dict:
    """BLACS-specific C migration (original behavior)."""
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

    # Patch Bdef.h with extended-precision type definitions and macros
    bdef_path = output_dir / 'Bdef.h'
    if bdef_path.exists():
        _patch_bdef_header(bdef_path, kind, template_vars)

    # Generate MPI requirement check for KIND=16
    if kind == 16:
        _generate_mpi_real16_check(output_dir)

    return {
        'cloned': cloned,
        'template_vars': template_vars,
    }


# C type underlying each target KIND.
_C_REAL_TYPE = {
    16: '__float128',
    10: 'long double',
}

# The 11 BLACS user-facing routine suffixes (same for each type prefix).
_BLACS_ROUTINE_SUFFIXES = [
    'gesd2d', 'gerv2d', 'gebs2d', 'gebr2d',
    'trsd2d', 'trrv2d', 'trbs2d', 'trbr2d',
    'gsum2d', 'gamx2d', 'gamn2d',
]


def _patch_bdef_header(bdef_path: Path, kind: int,
                       template_vars: dict[str, str]) -> None:
    """Add extended-precision type definitions and macros to Bdef.h."""
    rp = template_vars['RP']
    cp = template_vars['CP']
    real_type = template_vars['REAL_TYPE']      # e.g. QREAL
    complex_type = template_vars['COMPLEX_TYPE']  # e.g. XCOMPLEX
    c_type = _C_REAL_TYPE[kind]

    # --- Type definitions and prototypes ---
    type_block = f"""
/*
 *  Extended-precision types for migrated {{prefix}} routines.
 */
typedef {c_type} {real_type};
typedef struct {{{real_type} r, i;}} {complex_type};

void BI_{rp}mvcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
void BI_{rp}vmcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
""".replace('{prefix}', f'{rp}/{cp}')

    # --- Complex copy macros ---
    macro_block = f"""#define BI_{cp}mvcopy(m, n, A, lda, buff) \\
        BI_{rp}mvcopy(2*(m), (n), ({real_type} *) (A), 2*(lda), ({real_type} *) (buff))
#define BI_{cp}vmcopy(m, n, A, lda, buff) \\
        BI_{rp}vmcopy(2*(m), (n), ({real_type} *) (A), 2*(lda), ({real_type} *) (buff))
"""

    # --- Fortran name mangling defines ---
    def _mangling_block(prefix: str, transform) -> str:
        lines = []
        for suf in _BLACS_ROUTINE_SUFFIXES:
            src = f'{prefix}{suf}_'
            dst = transform(f'{prefix}{suf}')
            lines.append(f'#define {src:19s}{dst}')
        return '\n'.join(lines) + '\n'

    nochange_block = (_mangling_block(rp, lambda s: s) +
                      _mangling_block(cp, lambda s: s))
    upcase_block = (_mangling_block(rp, str.upper) +
                    _mangling_block(cp, str.upper))

    text = bdef_path.read_text(errors='replace')

    # Insert type definitions after the DCOMPLEX/SCOMPLEX typedefs
    text = text.replace(
        'typedef struct {float r, i;} SCOMPLEX;\n',
        'typedef struct {float r, i;} SCOMPLEX;\n' + type_block
    )

    # Insert copy macros after the existing BI_zvmcopy macro
    zvmcopy_line = '#define BI_zvmcopy(m, n, A, lda, buff) \\\n        BI_dvmcopy(2*(m), (n), (double *) (A), 2*(lda), (double *) (buff))'
    text = text.replace(zvmcopy_line, zvmcopy_line + '\n' + macro_block)

    # Insert name mangling defines in NOCHANGE section (after zgamn2d)
    nochange_marker = f'#define zgamn2d_   zgamn2d\n'
    text = text.replace(nochange_marker, nochange_marker + nochange_block)

    # Insert name mangling defines in UPCASE section (after ZGAMN2D)
    upcase_marker = f'#define zgamn2d_   ZGAMN2D\n'
    text = text.replace(upcase_marker, upcase_marker + upcase_block)

    bdef_path.write_text(text)


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
