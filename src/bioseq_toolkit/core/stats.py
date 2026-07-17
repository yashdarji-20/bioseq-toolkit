"""Summary statistics for DNA sequences.

Bundles the composition metrics plus molecular weight and melting temperature
into a single :class:`SequenceStats` result. All functions reuse the tested
``composition`` helpers and raise :class:`InvalidSequenceError` on non-DNA
input, matching the rest of the toolkit.
"""

from __future__ import annotations

from dataclasses import dataclass

from bioseq_toolkit.core.composition import (
    _require_valid_dna,
    at_content,
    gc_content,
    nucleotide_frequency,
)

# Average molecular weights (g/mol) of the four DNA deoxyribonucleotide
# monophosphates. Source: standard biochemistry reference values.
_DNA_NUCLEOTIDE_WEIGHTS: dict[str, float] = {
    "A": 331.2222,
    "C": 307.1971,
    "G": 347.2212,
    "T": 322.2085,
}

# Mass of water (g/mol) released for each phosphodiester bond formed when
# nucleotides are joined into a strand.
_WATER_WEIGHT = 18.0153

# The Wallace rule is only appropriate for short oligonucleotides.
_WALLACE_MAX_LENGTH = 14


@dataclass(frozen=True)
class SequenceStats:
    """A bundle of summary statistics for a DNA sequence.

    Attributes:
        length: Number of bases.
        gc_percent: GC content as a percentage.
        at_percent: AT content as a percentage.
        base_counts: Count of each of A, C, G, T.
    """

    length: int
    gc_percent: float
    at_percent: float
    base_counts: dict[str, int]


def molecular_weight(sequence: str) -> float:
    """Estimate the molecular weight of a single-stranded DNA sequence.

    Sums the weight of each nucleotide, then subtracts one water molecule per
    phosphodiester bond (``length - 1`` bonds) formed when the strand is built.

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        The estimated molecular weight in g/mol (Daltons).

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    total = sum(_DNA_NUCLEOTIDE_WEIGHTS[base] for base in seq)
    total -= _WATER_WEIGHT * (len(seq) - 1)
    return total


def melting_temperature(sequence: str) -> float:
    """Estimate melting temperature (Tm) in degrees Celsius (Wallace rule).

    Uses the Wallace rule ``Tm = 2*(A+T) + 4*(G+C)``, which is appropriate for
    short primers. A note is raised via exception only for invalid sequences;
    long sequences still compute but the estimate is less accurate.

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        The estimated melting temperature in degrees Celsius.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    a_t = seq.count("A") + seq.count("T")
    g_c = seq.count("G") + seq.count("C")
    return 2 * a_t + 4 * g_c


def sequence_stats(sequence: str) -> SequenceStats:
    """Compute a bundle of summary statistics for a DNA sequence.

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        A :class:`SequenceStats` with length, GC%, AT%, and base counts.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    return SequenceStats(
        length=len(seq),
        gc_percent=gc_content(seq),
        at_percent=at_content(seq),
        base_counts=nucleotide_frequency(seq),
    )
