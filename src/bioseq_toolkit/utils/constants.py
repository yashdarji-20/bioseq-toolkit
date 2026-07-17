"""Domain constants for BioSeq Toolkit.

This module is the single source of truth for biological facts used across
the toolkit: the kinds of sequences we support, their valid characters, and
the genetic code (codon -> amino acid). Every other module imports from here
rather than hardcoding letters like "ATCG", so a change is made in one place.
"""

from __future__ import annotations

from enum import Enum

# ---------------------------------------------------------------------------
# Sequence types
# ---------------------------------------------------------------------------


class SequenceType(str, Enum):
    """The kinds of biological sequence this toolkit understands.

    Subclassing ``str`` lets a member compare equal to its string value
    (``SequenceType.DNA == "DNA"``) and print nicely, while still giving us
    the safety and autocompletion of an enum.
    """

    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "PROTEIN"


# ---------------------------------------------------------------------------
# Valid characters for each sequence type
# ---------------------------------------------------------------------------
# frozenset -> immutable + fast membership testing ("A" in DNA_BASES).

DNA_BASES: frozenset[str] = frozenset("ACGT")
RNA_BASES: frozenset[str] = frozenset("ACGU")

# The 20 standard amino acids in single-letter form.
AMINO_ACIDS: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Convenience mapping so a validator can look up the right set by type.
VALID_CHARACTERS: dict[SequenceType, frozenset[str]] = {
    SequenceType.DNA: DNA_BASES,
    SequenceType.RNA: RNA_BASES,
    SequenceType.PROTEIN: AMINO_ACIDS,
}

# ---------------------------------------------------------------------------
# Complement rules (used later by reverse-complement)
# ---------------------------------------------------------------------------

DNA_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}
RNA_COMPLEMENT: dict[str, str] = {"A": "U", "U": "A", "C": "G", "G": "C"}

# ---------------------------------------------------------------------------
# Special codons
# ---------------------------------------------------------------------------

START_CODON: str = "AUG"
STOP_CODONS: frozenset[str] = frozenset({"UAA", "UAG", "UGA"})

# ---------------------------------------------------------------------------
# The genetic code: RNA codon -> amino acid single letter ("*" = stop)
# ---------------------------------------------------------------------------
# Written against RNA (U, not T) because translation acts on mRNA.

CODON_TABLE: dict[str, str] = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}