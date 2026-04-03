"""Precision prefix classification and rename map generation.

Supports two prefix styles:
  - "direct":    S/D/C/Z is the first character (BLAS, LAPACK)
  - "scalapack": P + S/D/C/Z (PDGESV, PZHEEV)
"""

# Target prefix mapping by KIND
PREFIX_MAP: dict[int, dict[str, str]] = {
    10: {'S': 'E', 'D': 'E', 'C': 'Y', 'Z': 'Y'},
    16: {'S': 'Q', 'D': 'Q', 'C': 'X', 'Z': 'X'},
}

PRECISION_CHARS = frozenset('SDCZ')


def classify_prefix(name: str, style: str = 'direct') -> tuple[str, str]:
    """Classify a symbol's precision prefix.

    Returns (prefix_char, base_name) where prefix_char is one of
    S, D, C, Z, or '' for non-precision symbols.

    Args:
        name: Uppercase symbol name (e.g., 'DGEMM', 'PDGESV', 'LSAME')
        style: 'direct' for BLAS/LAPACK, 'scalapack' for P-prefixed
    """
    if len(name) < 2:
        return ('', name)

    if style == 'scalapack':
        if len(name) > 2 and name[0] == 'P' and name[1] in PRECISION_CHARS:
            return (name[1], 'P' + name[2:])
        # Fall through to direct classification for non-P-prefixed
        # routines in ScaLAPACK (e.g., helper routines)

    if name[0] in PRECISION_CHARS:
        return (name[0], name[1:])

    return ('', name)


def build_rename_map(symbols: set[str], target_kind: int,
                     style: str = 'direct') -> dict[str, str]:
    """Build old_name → new_name mapping.

    Only renames symbols that have at least one precision-prefixed variant
    in the symbol set (avoids renaming LSAME, XERBLA, etc.).
    """
    pmap = PREFIX_MAP.get(target_kind)
    if pmap is None:
        raise ValueError(f'Unsupported target kind: {target_kind}')

    # Collect base names that have at least one precision variant
    bases_with_variants: set[str] = set()
    for sym in symbols:
        pfx, base = classify_prefix(sym, style)
        if pfx:
            bases_with_variants.add(base)

    rename: dict[str, str] = {}
    for sym in symbols:
        pfx, base = classify_prefix(sym, style)
        if pfx and base in bases_with_variants:
            new_pfx = pmap[pfx]
            if style == 'scalapack' and base.startswith('P'):
                # P + new_prefix + rest
                rename[sym] = 'P' + new_pfx + base[1:]
            else:
                rename[sym] = new_pfx + base
    return rename


def target_prefix(kind: int, is_complex: bool) -> str:
    """Return the target prefix character for a given KIND and type."""
    pmap = PREFIX_MAP[kind]
    if is_complex:
        return pmap['Z']  # Z and C both map to the same target
    else:
        return pmap['D']  # D and S both map to the same target
