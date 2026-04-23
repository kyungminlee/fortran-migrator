"""Identify DOUBLE PRECISION declarations that should be preserved
(not promoted) during MUMPS type migration.

Rule: for each DOUBLE PRECISION declaration in a d*/z* file, check whether
the same declaration (modulo whitespace) appears in the paired s*/c* file.
If yes, it is intentionally DP on both sides (timing, flop counters, memory
estimates, or anything interfacing with MPI_WTIME) and must be left alone;
promoting it would create divergence with the s-migrated output and would
break the MPI_WTIME return-kind contract.

Usage:
    python mumps_find_keep_kind.py <s_file> <d_file>

Output (stdout):
    <d_file>:<line_number>    <original line>
"""

import re
import sys
from pathlib import Path

DP_RE = re.compile(r"^\s*DOUBLE\s+PRECISION\b", re.IGNORECASE)


def normalize(line: str) -> str:
    """Collapse whitespace and lowercase for comparison."""
    # Strip trailing newline and any column-73+ sequence numbers (fixed-form).
    body = line.rstrip("\n")[:72]
    return re.sub(r"\s+", " ", body).strip().lower()


def dp_lines(path: Path) -> dict[str, list[int]]:
    """Return {normalized_line: [line_numbers]} for every DP declaration."""
    out: dict[str, list[int]] = {}
    for i, line in enumerate(path.read_text().splitlines(), start=1):
        if DP_RE.match(line):
            out.setdefault(normalize(line), []).append(i)
    return out


def main(s_path: Path, d_path: Path) -> int:
    s_set = set(dp_lines(s_path).keys())
    d_map = dp_lines(d_path)
    d_text = d_path.read_text().splitlines()

    keep = []
    promote = []
    for norm, lines in d_map.items():
        for ln in lines:
            entry = (ln, d_text[ln - 1])
            (keep if norm in s_set else promote).append(entry)

    keep.sort()
    promote.sort()

    print(f"# {d_path}")
    print(f"# keep-kind: {len(keep)}   promote: {len(promote)}")
    print("# --- keep-kind (DP in both s and d — preserve) ---")
    for ln, text in keep:
        print(f"{d_path}:{ln}\t{text.rstrip()}")
    print("# --- promote (DP only in d — working precision) ---")
    for ln, text in promote:
        print(f"{d_path}:{ln}\t{text.rstrip()}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
