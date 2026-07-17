"""Translation (DNA/RNA -> protein) for BioSeq Toolkit.

Translation reads a sequence in codons (groups of three bases) and maps each
codon to an amino acid via :data:`CODON_TABLE`. This module composes earlier
building blocks: it validates input, transcribes DNA to RNA, then translates.

Design choices (documented so behaviour is predictable):
    * Input is treated as DNA and transcribed to RNA first, because the codon
      table is keyed by RNA codons.
    * Translation starts at position 0 (no start-codon search).
    * A trailing partial codon (1-2 leftover bases) is ignored.
    * ``to_stop=True`` (default) ends translation at the first stop codon;
      ``to_stop=False`` translates the whole sequence, marking stops as ``*``.
"""

from __future__ import annotations

from bioseq_toolkit.core.composition import _require_valid_dna
from bioseq_toolkit.core.transforms import transcribe
from bioseq_toolkit.utils.constants import CODON_TABLE

_CODON_LENGTH = 3


def translate(sequence: str, *, to_stop: bool = True) -> str:
    """Translate a DNA sequence into a protein (amino acid string).

    The sequence is validated, transcribed to RNA, then read codon by codon.

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).
        to_stop: If ``True`` (default), stop translating at the first stop
            codon and do not include it. If ``False``, translate the entire
            sequence, representing each stop codon as ``"*"``.

    Returns:
        The protein as a string of single-letter amino acid codes. May be
        empty if translation stops immediately at a leading stop codon.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    _require_valid_dna(sequence)  # raises on invalid input
    rna = transcribe(sequence)

    amino_acids: list[str] = []
    # Subtracting (_CODON_LENGTH - 1) makes the loop stop before any partial
    # trailing codon, so leftover 1-2 bases are ignored automatically.
    for i in range(0, len(rna) - (_CODON_LENGTH - 1), _CODON_LENGTH):
        codon = rna[i : i + _CODON_LENGTH]
        amino_acid = CODON_TABLE[codon]
        if amino_acid == "*" and to_stop:
            break
        amino_acids.append(amino_acid)

    return "".join(amino_acids)
