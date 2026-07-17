"""Nucleotide composition analysis for BioSeq Toolkit.

Provides GC content, AT content, and per-base frequency for DNA sequences.
These functions RAISE :class:`InvalidSequenceError` on non-DNA input rather
than returning a misleading number, because there is no meaningful GC content
for a sequence that is not DNA.
"""

from __future__ import annotations

from collections import Counter

from bioseq_toolkit.core.validation import _normalize, is_valid_dna
from bioseq_toolkit.utils.exceptions import InvalidSequenceError


def _require_valid_dna(sequence: str) -> str:
    """Normalize and validate a DNA sequence, or raise.

    Args:
        sequence: The raw input sequence.

    Returns:
        The normalized (stripped, uppercased) sequence.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    if not is_valid_dna(sequence):
        raise InvalidSequenceError(f"Not a valid DNA sequence: {sequence!r}")
    return _normalize(sequence)


def gc_content(sequence: str) -> float:
    """Return the GC content of a DNA sequence as a percentage (0-100).

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        The percentage of bases that are G or C, as a float.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100


def at_content(sequence: str) -> float:
    """Return the AT content of a DNA sequence as a percentage (0-100).

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        The percentage of bases that are A or T, as a float.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    at = seq.count("A") + seq.count("T")
    return at / len(seq) * 100


def nucleotide_frequency(sequence: str) -> dict[str, int]:
    """Return a count of each base in a DNA sequence.

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        A dict mapping each of A, C, G, T to its count. Bases that do not
        appear are included with a count of 0, so the shape is predictable.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    counts = Counter(seq)
    # Guarantee all four keys exist, even when a base is absent.
    return {base: counts.get(base, 0) for base in "ACGT"}