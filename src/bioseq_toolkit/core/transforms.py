"""Sequence transformations for BioSeq Toolkit.

Implements the two most fundamental DNA transformations:

* :func:`reverse_complement` - the opposite strand read 5'->3'.
* :func:`transcribe` - DNA to messenger RNA (every ``T`` becomes ``U``).

Both raise :class:`InvalidSequenceError` on non-DNA input, matching the
contract used elsewhere in the toolkit.
"""

from __future__ import annotations

from bioseq_toolkit.core.composition import _require_valid_dna

# Translation tables are built once at import time (not per call) for speed.
# maketrans("ACGT", "TGCA") maps A<->T and C<->G in a single C-level pass.
_COMPLEMENT_TABLE = str.maketrans("ACGT", "TGCA")
_TRANSCRIBE_TABLE = str.maketrans("T", "U")


def complement(sequence: str) -> str:
    """Return the base-by-base complement of a DNA sequence (not reversed).

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        The complement, e.g. ``"ATCG"`` -> ``"TAGC"``.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    return seq.translate(_COMPLEMENT_TABLE)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    The reverse complement is the opposite strand written in the biological
    5'->3' reading direction: complement every base, then reverse the string.

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        The reverse complement, e.g. ``"ATCG"`` -> ``"CGAT"``.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    # [::-1] reverses the string; .translate applies the complement map.
    return seq.translate(_COMPLEMENT_TABLE)[::-1]


def transcribe(sequence: str) -> str:
    """Transcribe DNA into messenger RNA (replace every ``T`` with ``U``).

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        The RNA sequence, e.g. ``"ATCG"`` -> ``"AUCG"``.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    return seq.translate(_TRANSCRIBE_TABLE)