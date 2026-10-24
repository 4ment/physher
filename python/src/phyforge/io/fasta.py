"""Minimal FASTA reader.

Deliberately dependency-free: phyforge only ever needs taxon names and, when
embedding an alignment, the raw sequence strings.
"""

from __future__ import annotations

from typing import List, NamedTuple


class Sequence(NamedTuple):
    taxon: str
    sequence: str


def read_fasta(path: str) -> List[Sequence]:
    """Read ``path``, preserving the order the sequences appear in."""
    sequences: List[Sequence] = []
    taxon = None
    chunks: List[str] = []
    with open(path, "r") as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if taxon is not None:
                    sequences.append(Sequence(taxon, "".join(chunks)))
                # The name is everything up to the first whitespace, which is
                # how physher's own reader splits it.
                taxon = line[1:].split()[0] if line[1:].split() else ""
                chunks = []
            else:
                chunks.append(line)
    if taxon is not None:
        sequences.append(Sequence(taxon, "".join(chunks)))
    if not sequences:
        raise ValueError(f"no sequences found in {path}")
    return sequences
