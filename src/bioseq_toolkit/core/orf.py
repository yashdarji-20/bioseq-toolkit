"""Open Reading Frame (ORF) detection for BioSeq Toolkit.

An ORF is a stretch of DNA that could encode a protein: it begins with the
start codon ``ATG``, proceeds in codons, and ends at an in-frame stop codon
(``TAA``, ``TAG``, ``TGA``).

Because a gene may sit in any of three reading frames on either strand, this
module performs the standard *six-frame* search: frames 0-2 on the forward
strand and frames 0-2 on the reverse complement.
"""

from __future__ import annotations

from dataclasses import dataclass

from bioseq_toolkit.core.composition import _require_valid_dna
from bioseq_toolkit.core.transforms import reverse_complement
from bioseq_toolkit.core.translation import translate

_CODON_LENGTH = 3
_START_CODON_DNA = "ATG"
_STOP_CODONS_DNA = frozenset({"TAA", "TAG", "TGA"})

# Minimum ORF length in bases. Short ORFs arise by chance very often, so a
# floor keeps results meaningful. 75 bases (25 codons) is a common default.
_DEFAULT_MIN_LENGTH = 75


@dataclass(frozen=True)
class ORF:
    """A single open reading frame.

    Attributes:
        strand: ``"+"`` for the forward strand, ``"-"`` for the reverse.
        frame: Reading frame offset (0, 1, or 2) within that strand.
        start: 0-based start index of the ORF, in the coordinates of the
            strand it was found on.
        end: 0-based exclusive end index (just past the stop codon).
        sequence: The ORF's DNA sequence, including start and stop codons.
        protein: The translated protein (stop codon not included).
    """

    strand: str
    frame: int
    start: int
    end: int
    sequence: str
    protein: str

    @property
    def length(self) -> int:
        """Length of the ORF in bases."""
        return len(self.sequence)


def _find_orfs_in_strand(seq: str, strand: str, min_length: int) -> list[ORF]:
    """Find ORFs in all three reading frames of one strand.

    Args:
        seq: The (already validated and normalized) strand sequence.
        strand: ``"+"`` or ``"-"``, recorded on each result.
        min_length: Minimum ORF length in bases.

    Returns:
        A list of :class:`ORF` found on this strand.
    """
    orfs: list[ORF] = []

    for frame in range(_CODON_LENGTH):
        # Walk this frame codon by codon looking for start codons.
        for i in range(frame, len(seq) - (_CODON_LENGTH - 1), _CODON_LENGTH):
            if seq[i : i + _CODON_LENGTH] != _START_CODON_DNA:
                continue

            # Found a start: scan forward in-frame for a stop codon.
            for j in range(i, len(seq) - (_CODON_LENGTH - 1), _CODON_LENGTH):
                codon = seq[j : j + _CODON_LENGTH]
                if codon in _STOP_CODONS_DNA:
                    end = j + _CODON_LENGTH  # exclusive, includes the stop
                    orf_seq = seq[i:end]
                    if len(orf_seq) >= min_length:
                        orfs.append(
                            ORF(
                                strand=strand,
                                frame=frame,
                                start=i,
                                end=end,
                                sequence=orf_seq,
                                protein=translate(orf_seq),
                            )
                        )
                    break  # this start is resolved; move to the next one

    return orfs


def find_orfs(
    sequence: str,
    *,
    min_length: int = _DEFAULT_MIN_LENGTH,
    both_strands: bool = True,
) -> list[ORF]:
    """Find open reading frames in a DNA sequence.

    Searches all three reading frames on the forward strand and, by default,
    all three on the reverse complement (a six-frame search).

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).
        min_length: Minimum ORF length in bases; shorter ORFs are discarded.
        both_strands: If ``True`` (default), also search the reverse strand.

    Returns:
        A list of :class:`ORF`, sorted longest first.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)

    orfs = _find_orfs_in_strand(seq, "+", min_length)
    if both_strands:
        orfs += _find_orfs_in_strand(reverse_complement(seq), "-", min_length)

    # Longest first makes the most promising candidates easy to spot.
    orfs.sort(key=lambda orf: orf.length, reverse=True)
    return orfs
