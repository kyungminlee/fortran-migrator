"""Precision prefix classification via empirical pattern discovery.

Discovers precision structure by analyzing the full symbol set in two passes:

  Pass 1 (single precision): Replace S→@R@, C→@C@ at candidate positions.
  Pass 2 (double precision): Replace D→@R@, Z→@C@ at candidate positions.

Within each pass, symbols are grouped by their normalized pattern.
Groups of 2+ reveal precision-dependent families. Families from both
passes that share the same base pattern are merged into a single family.

This naturally handles any prefix convention without hardcoding:
  - Direct:      SGEMM/DGEMM/CGEMM/ZGEMM → @T@GEMM
  - I-prefix:    ISAMAX/IDAMAX/ICAMAX/IZAMAX → I@T@AMAX
  - ScaLAPACK:   PSGESV/PDGESV → P@T@GESV
  - Cross-type:  CSROT(pass1)→@@ROT, ZDROT(pass2)→@@ROT → merged
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

# Target prefix mapping by KIND
PREFIX_MAP: dict[int, dict[str, str]] = {
    10: {'S': 'E', 'D': 'E', 'C': 'Y', 'Z': 'Y'},
    16: {'S': 'Q', 'D': 'Q', 'C': 'X', 'Z': 'X'},
}

# Two passes: single-precision chars and double-precision chars
SINGLE_CHARS = {'S', 'C'}  # S=real, C=complex
DOUBLE_CHARS = {'D', 'Z'}  # D=real, Z=complex


@dataclass
class PrecisionFamily:
    """A group of symbols sharing the same base pattern with
    precision-character variations."""
    pattern: str                    # e.g., '@GEMM', 'I@AMAX', '@@ROT'
    members: dict[str, str]         # precision_key → symbol name
    prefix_positions: list[int]     # positions of precision chars


@dataclass
class SymbolClassification:
    """Result of analyzing a set of symbols for precision patterns."""
    families: list[PrecisionFamily] = field(default_factory=list)
    independent: set[str] = field(default_factory=set)
    _symbol_to_family: dict[str, PrecisionFamily] = field(
        default_factory=dict, repr=False)

    def is_precision_dependent(self, name: str) -> bool:
        return name.upper() in self._symbol_to_family

    def get_family(self, name: str) -> PrecisionFamily | None:
        return self._symbol_to_family.get(name.upper())

    def target_name(self, name: str, target_kind: int) -> str | None:
        """Compute the target name for a symbol at the given KIND."""
        upper = name.upper()
        family = self._symbol_to_family.get(upper)
        if family is None:
            return None

        pmap = PREFIX_MAP[target_kind]
        result = list(upper)
        for pos in family.prefix_positions:
            old_char = upper[pos]
            if old_char in pmap:
                result[pos] = pmap[old_char]
        return ''.join(result)

    def build_rename_map(self, target_kind: int) -> dict[str, str]:
        """Build old_name → new_name mapping for all precision-dependent symbols."""
        rename: dict[str, str] = {}
        for sym in self._symbol_to_family:
            new_name = self.target_name(sym, target_kind)
            if new_name and new_name != sym:
                rename[sym] = new_name
        return rename


def _find_families_single_pass(
    symbols: set[str],
    replacement_chars: set[str],
) -> dict[tuple[str, tuple[int, ...]], set[str]]:
    """Run one pass of pattern discovery.

    For each symbol, find positions containing chars from replacement_chars.
    Try all non-empty subsets of those positions as replacement candidates.
    Group symbols by their normalized pattern.

    Returns: { (pattern, positions) → set of matching symbols }
    """
    groups: dict[tuple[str, tuple[int, ...]], set[str]] = {}

    for sym in symbols:
        # Find candidate positions
        candidate_positions = [
            i for i, ch in enumerate(sym) if ch in replacement_chars
        ]
        if not candidate_positions:
            continue

        # Try all non-empty subsets of candidate positions (limit to 3)
        for r in range(1, min(len(candidate_positions) + 1, 4)):
            for positions in combinations(candidate_positions, r):
                # Build pattern by replacing selected positions with '@'
                chars = list(sym)
                for pos in positions:
                    chars[pos] = '@'
                pattern = ''.join(chars)

                key = (pattern, positions)
                if key not in groups:
                    groups[key] = set()
                groups[key].add(sym)

    return groups


def classify_symbols(symbols: set[str]) -> SymbolClassification:
    """Analyze a set of symbols to discover precision families.

    Two-pass approach:
      Pass 1: Replace S and C at candidate positions → find single-precision families
      Pass 2: Replace D and Z at candidate positions → find double-precision families
      Merge: Families from both passes with the same base pattern are merged.
    """
    upper_symbols = {s.upper() for s in symbols}

    # Run both passes
    single_groups = _find_families_single_pass(upper_symbols, SINGLE_CHARS)
    double_groups = _find_families_single_pass(upper_symbols, DOUBLE_CHARS)

    # Merge: groups from both passes with the same (pattern, positions)
    # key are combined. A family is valid if it has 2+ members from
    # EITHER pass or from the merged result.
    all_groups: dict[tuple[str, tuple[int, ...]], set[str]] = {}

    for key, members in single_groups.items():
        all_groups[key] = members

    for key, members in double_groups.items():
        if key in all_groups:
            all_groups[key] = all_groups[key] | members
        else:
            all_groups[key] = members

    # Keep groups with 2+ members
    candidates: dict[tuple[str, tuple[int, ...]], set[str]] = {
        key: members for key, members in all_groups.items()
        if len(members) >= 2
    }

    # Sort candidates: most members first, then fewest positions (most specific)
    sorted_candidates = sorted(
        candidates.items(),
        key=lambda item: (-len(item[1]), len(item[0][1]))
    )

    # Greedily assign symbols to families
    assigned: set[str] = set()
    result = SymbolClassification()

    for (pattern, positions), members in sorted_candidates:
        unassigned = members - assigned
        if len(unassigned) < 2:
            continue

        # Build member dict: precision_key → symbol
        member_dict: dict[str, str] = {}
        for sym in sorted(unassigned):
            prec_key = ''.join(sym[p] for p in positions)
            member_dict[prec_key] = sym

        family = PrecisionFamily(
            pattern=pattern,
            members=member_dict,
            prefix_positions=list(positions),
        )
        result.families.append(family)

        for sym in unassigned:
            assigned.add(sym)
            result._symbol_to_family[sym] = family

    # Everything not assigned is independent
    result.independent = upper_symbols - assigned
    return result


# Convenience wrappers

def build_rename_map(symbols: set[str], target_kind: int,
                     style: str = 'direct') -> dict[str, str]:
    """Build old_name → new_name mapping.

    The style parameter is accepted for backward compatibility but
    ignored — the classifier discovers the prefix structure empirically.
    """
    classification = classify_symbols(symbols)
    return classification.build_rename_map(target_kind)


def target_prefix(kind: int, is_complex: bool) -> str:
    """Return the target prefix character for a given KIND and type."""
    pmap = PREFIX_MAP[kind]
    return pmap['Z'] if is_complex else pmap['D']
