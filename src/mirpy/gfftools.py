from __future__ import annotations
from pathlib import Path
import gzip
from typing import Iterable, TextIO, Collection, Union


def _open_text_auto(path: str | Path, mode: str = "rt") -> TextIO:
    p = Path(path)
    if p.suffix.lower() == ".gz":
        return gzip.open(p, mode, encoding="utf-8", errors="replace")
    return open(p, mode, encoding="utf-8", errors="replace")


def subset_gff_by_criteria(
        in_path: str | Path,
        out_path: str | Path,
        col_index: int = 2,
        split_char: str = "\t",
        criteria: Union[str, Collection[str]] = "miRNA",
        case_insensitive: bool = False,
) -> int:
    """
    Copy GFF3 header lines and only those rows whose Nth column equals criteria.
    Works with .gff3 or .gff3.gz for input and output (based on extension).
    """

    # Process criteria
    if isinstance(criteria, str):
        crit_set = {criteria}
    else:
        crit_set = set(criteria)

    if case_insensitive:
        crit_set = {c.lower() for c in crit_set}

    # Input GFF3
    fin = _open_text_auto(in_path, "rt")
    # Output GFF3
    pout = Path(out_path)
    if pout.suffix.lower() == ".gz":
        fout = gzip.open(pout, "wt", encoding="utf-8")
    else:
        fout = open(pout, "wt", encoding="utf-8")

    # Number of output rows
    count = 0
    # Iterate the input and write out to output
    with fin, fout:
        for raw in fin:
            if not raw.strip():
                continue
            if raw.startswith("#"):
                fout.write(raw)
                continue

            cols = raw.rstrip("\n").split(split_char)
            if col_index < 0 or col_index >= len(cols):
                # malformed line; skip safely
                continue

            cell = cols[col_index]
            key = cell.lower() if case_insensitive else cell

            if key in crit_set:
                fout.write(raw)
                count += 1

    # Return written line count
    return count