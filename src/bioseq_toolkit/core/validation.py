"""Sequence validation for BioSeq Toolkit.

Answers the question "is this a well-formed DNA / RNA / protein sequence?"
All three validators share one private core (:func:`_is_valid`) so the rule
"every character must be in the allowed alphabet" lives in exactly one place.

Validation rules (Version 1):
    * Surrounding whitespace is ignored and input is treated case-insensitively.
    * An empty sequence (after stripping) is NOT valid.
    * Only canonical characters are allowed: A/C/G/T for DNA, A/C/G/U for RNA,
      the 20 standard amino acids for protein. Ambiguity codes (e.g. ``N``)
      are rejected in this version.
"""

from __future__ import annotations

from bioseq_toolkit.utils.constants import (
    AMINO_ACIDS,
    DNA_BASES,
    RNA_BASES,
    SequenceType,
    VALID_CHARACTERS,
)


def _normalize(sequence: str) -> str:
    """Strip surrounding whitespace and uppercase a sequence.

    Args:
        sequence: The raw input sequence.

    Returns:
        The cleaned sequence, ready for validation.
    """
    return sequence.strip().upper()


def _is_valid(sequence: str, sequence_type: SequenceType) -> bool:
    """Core validation shared by all sequence types.

    Args:
        sequence: The raw input sequence.
        sequence_type: Which alphabet to validate against.

    Returns:
        ``True`` if the normalized sequence is non-empty and every character
        belongs to the alphabet for ``sequence_type``; otherwise ``False``.
    """
    normalized = _normalize(sequence)
    if not normalized:
        return False
    allowed = VALID_CHARACTERS[sequence_type]
    # set(normalized) is the unique characters used; if that set has nothing
    # outside `allowed`, the sequence is valid. issubset expresses this directly.
    return set(normalized).issubset(allowed)


def is_valid_dna(sequence: str) -> bool:
    """Return ``True`` if ``sequence`` is valid DNA (only A/C/G/T)."""
    return _is_valid(sequence, SequenceType.DNA)


def is_valid_rna(sequence: str) -> bool:
    """Return ``True`` if ``sequence`` is valid RNA (only A/C/G/U)."""
    return _is_valid(sequence, SequenceType.RNA)


def is_valid_protein(sequence: str) -> bool:
    """Return ``True`` if ``sequence`` is a valid protein (20 amino acids)."""
    return _is_valid(sequence, SequenceType.PROTEIN)