from dataclasses import dataclass

# A more memory efficient way of storing just essential information from all the alignment info
@dataclass
class AlignmentData:
    """Minimal alignment data to reduce memory footprint."""
    __slots__ = ('chr', 'start_1b', 'end_1b', 'is_reverse', 'nm')
    chr: str
    start_1b: int
    end_1b: int
    is_reverse: bool
    nm: int

# Mature miRNA
@dataclass
class Mature:
    name: str
    chr: str
    strand: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
