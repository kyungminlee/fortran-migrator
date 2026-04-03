"""Precision prefix classification via empirical pattern discovery.

Discovers precision structure by analyzing the full symbol set in two passes:

  Pass 1 (single precision): Replace S→@, C→@ at candidate positions.
  Pass 2 (double precision): Replace D→@, Z→@ at candidate positions.

Within each pass, symbols are grouped by their normalized pattern.
Groups of 2+ reveal precision-dependent families. Families from both
passes sharing the same base pattern are merged.

Each replaced position carries a **type tag** (R=real, C=complex) so
that multi-position patterns like DZNRM2→@@NRM2 know position 0 is
real and position 1 is complex — not just "some precision char."

This naturally handles any prefix convention without hardcoding:
  - Direct:      SGEMM/DGEMM/CGEMM/ZGEMM → @GEMM
  - I-prefix:    ISAMAX/IDAMAX/ICAMAX/IZAMAX → I@AMAX
  - ScaLAPACK:   PSGESV/PDGESV → P@GESV
  - Cross-type:  CSROT/ZDROT → @@ROT (pos0=C, pos1=R)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

# Target prefix mapping by KIND
PREFIX_MAP: dict[int, dict[str, str]] = {
    10: {'S': 'E', 'D': 'E', 'C': 'Y', 'Z': 'Y'},
    16: {'S': 'Q', 'D': 'Q', 'C': 'X', 'Z': 'X'},
}

# Character → type tag
CHAR_TYPE: dict[str, str] = {
    'S': 'R', 'D': 'R',    # real
    'C': 'C', 'Z': 'C',    # complex
}

# Two passes: single-precision chars and double-precision chars
SINGLE_CHARS = {'S', 'C'}
DOUBLE_CHARS = {'D', 'Z'}


@dataclass
class PrecisionSlot:
    """One position in a symbol that carries a precision character."""
    position: int
    type_tag: str    # 'R' (real) or 'C' (complex)


@dataclass
class PrecisionFamily:
    """A group of symbols sharing the same base pattern with
    precision-character variations."""
    pattern: str                     # e.g., '@GEMM', '@@ROT'
    members: dict[str, str]          # precision_key → symbol name
    slots: list[PrecisionSlot]       # positions + type tags

    @property
    def prefix_positions(self) -> list[int]:
        return [s.position for s in self.slots]


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
        for slot in family.slots:
            old_char = upper[slot.position]
            if old_char in pmap:
                result[slot.position] = pmap[old_char]
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
) -> dict[tuple[str, tuple[int, ...]], dict[str, tuple[str, ...]]]:
    """Run one pass of pattern discovery.

    Returns: { (pattern, positions) → { symbol: type_tags_tuple } }
    where type_tags_tuple is e.g. ('R',) or ('C', 'R') for each position.
    """
    groups: dict[tuple[str, tuple[int, ...]], dict[str, tuple[str, ...]]] = {}

    for sym in symbols:
        # Find candidate positions
        candidate_positions = [
            i for i, ch in enumerate(sym) if ch in replacement_chars
        ]
        if not candidate_positions:
            continue

        # Try all non-empty subsets (limit to 3 positions)
        for r in range(1, min(len(candidate_positions) + 1, 4)):
            for positions in combinations(candidate_positions, r):
                # Build pattern
                chars = list(sym)
                type_tags = []
                for pos in positions:
                    type_tags.append(CHAR_TYPE[sym[pos]])
                    chars[pos] = '@'
                pattern = ''.join(chars)

                key = (pattern, positions)
                if key not in groups:
                    groups[key] = {}
                groups[key][sym] = tuple(type_tags)

    return groups


def classify_symbols(symbols: set[str]) -> SymbolClassification:
    """Analyze a set of symbols to discover precision families.

    Two-pass approach:
      Pass 1: Replace S and C → find single-precision families
      Pass 2: Replace D and Z → find double-precision families
      Merge: Families from both passes with the same pattern are merged,
             but only if their type tags are consistent (same R/C at
             each position).
    """
    upper_symbols = {s.upper() for s in symbols}

    # Run both passes
    single_groups = _find_families_single_pass(upper_symbols, SINGLE_CHARS)
    double_groups = _find_families_single_pass(upper_symbols, DOUBLE_CHARS)

    # Merge: combine groups from both passes with the same (pattern, positions)
    # Only merge if type tags at each position are consistent.
    all_groups: dict[tuple[str, tuple[int, ...]], dict[str, tuple[str, ...]]] = {}

    for key, members in single_groups.items():
        all_groups[key] = dict(members)

    for key, members in double_groups.items():
        if key in all_groups:
            # Merge passes. Type tags may differ at a position when
            # one pass sees R (S or D) and the other sees C (C or Z).
            # This is fine — it means the slot carries both real and
            # complex variants (a "T" slot, like in a SDCZ quartet).
            # We upgrade the tag to 'T' when merging R+C.
            existing = all_groups[key]
            existing_tags = next(iter(existing.values()))
            new_tags = next(iter(members.values()))

            # Compute merged type tags
            merged_tags = tuple(
                'T' if et != nt else et
                for et, nt in zip(existing_tags, new_tags)
            )
            # Update all existing members to use merged tags
            for sym in list(existing):
                existing[sym] = merged_tags
            # Add new members with merged tags
            for sym in members:
                existing[sym] = merged_tags
        else:
            all_groups[key] = dict(members)

    # Keep groups with 2+ members
    candidates = {
        key: members for key, members in all_groups.items()
        if len(members) >= 2
    }

    # Sort: most members first, then fewest positions
    sorted_candidates = sorted(
        candidates.items(),
        key=lambda item: (-len(item[1]), len(item[0][1]))
    )

    # Greedily assign symbols to families
    assigned: set[str] = set()
    result = SymbolClassification()

    for (pattern, positions), members in sorted_candidates:
        unassigned = {s: t for s, t in members.items() if s not in assigned}
        if len(unassigned) < 2:
            continue

        # Get type tags from any member (they're all consistent)
        type_tags = next(iter(unassigned.values()))

        # Build slots with type tags
        slots = [
            PrecisionSlot(position=pos, type_tag=tag)
            for pos, tag in zip(positions, type_tags)
        ]

        # Build member dict: precision_key → symbol
        member_dict: dict[str, str] = {}
        for sym in sorted(unassigned):
            prec_key = ''.join(sym[p] for p in positions)
            member_dict[prec_key] = sym

        family = PrecisionFamily(
            pattern=pattern,
            members=member_dict,
            slots=slots,
        )
        result.families.append(family)

        for sym in unassigned:
            assigned.add(sym)
            result._symbol_to_family[sym] = family

    result.independent = upper_symbols - assigned
    return result


# Convenience wrappers

def build_rename_map(symbols: set[str], target_kind: int,
                     style: str = 'direct') -> dict[str, str]:
    """Build old_name → new_name mapping."""
    classification = classify_symbols(symbols)
    return classification.build_rename_map(target_kind)


def target_prefix(kind: int, is_complex: bool) -> str:
    """Return the target prefix character for a given KIND and type."""
    pmap = PREFIX_MAP[kind]
    return pmap['Z'] if is_complex else pmap['D']
