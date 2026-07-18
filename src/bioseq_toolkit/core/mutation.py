"""SNP detection and mutation comparison for BioSeq Toolkit.

Compares two DNA sequences and reports where they differ. Sequences of equal
length are compared position by position; sequences of differing length are
first aligned with Needleman-Wunsch (see :mod:`bioseq_toolkit.core.alignment`)
so that insertions and deletions are handled correctly.

Substitutions are classified as *transitions* (purine<->purine or
pyrimidine<->pyrimidine, e.g. A->G) or *transversions* (purine<->pyrimidine,
e.g. A->C). Transitions occur roughly twice as often in nature, so the
transition/transversion (Ti/Tv) ratio is a standard quality metric in variant
analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from bioseq_toolkit.core.alignment import align
from bioseq_toolkit.core.composition import _require_valid_dna

_PURINES = frozenset("AG")
_PYRIMIDINES = frozenset("CT")
_GAP_CHARACTER = "-"


class MutationType(str, Enum):
    """The kind of difference observed at a position."""

    TRANSITION = "transition"
    TRANSVERSION = "transversion"
    INSERTION = "insertion"
    DELETION = "deletion"


@dataclass(frozen=True)
class Variant:
    """A single difference between two sequences.

    Attributes:
        position: 0-based position in the comparison (alignment coordinates
            when the sequences were aligned).
        reference: The base in the reference sequence, or ``"-"`` for an
            insertion.
        sample: The base in the sample sequence, or ``"-"`` for a deletion.
        mutation_type: How this difference is classified.
    """

    position: int
    reference: str
    sample: str
    mutation_type: MutationType

    def __str__(self) -> str:
        """Return a compact ``pos:REF>SAMPLE`` description."""
        return f"{self.position}:{self.reference}>{self.sample}"


@dataclass(frozen=True)
class MutationReport:
    """Summary of all differences between two sequences.

    Attributes:
        variants: Every difference found, in position order.
        aligned_reference: The reference sequence as compared (gapped if
            alignment was performed).
        aligned_sample: The sample sequence as compared.
    """

    variants: list[Variant]
    aligned_reference: str
    aligned_sample: str

    @property
    def total_variants(self) -> int:
        """Total number of differences found."""
        return len(self.variants)

    @property
    def transitions(self) -> int:
        """Number of transition substitutions."""
        return sum(
            1 for v in self.variants if v.mutation_type is MutationType.TRANSITION
        )

    @property
    def transversions(self) -> int:
        """Number of transversion substitutions."""
        return sum(
            1 for v in self.variants if v.mutation_type is MutationType.TRANSVERSION
        )

    @property
    def indels(self) -> int:
        """Number of insertions and deletions."""
        return sum(
            1
            for v in self.variants
            if v.mutation_type in (MutationType.INSERTION, MutationType.DELETION)
        )

    @property
    def ti_tv_ratio(self) -> float:
        """Transition/transversion ratio (0.0 when there are no transversions)."""
        if self.transversions == 0:
            return 0.0
        return self.transitions / self.transversions

    @property
    def identity(self) -> float:
        """Percentage of compared positions that are identical."""
        if not self.aligned_reference:
            return 0.0
        length = len(self.aligned_reference)
        return (length - self.total_variants) / length * 100


def _classify_substitution(reference: str, sample: str) -> MutationType:
    """Classify a substitution as a transition or a transversion.

    Args:
        reference: The reference base.
        sample: The sample base.

    Returns:
        :attr:`MutationType.TRANSITION` if both bases are purines or both are
        pyrimidines, otherwise :attr:`MutationType.TRANSVERSION`.
    """
    both_purines = reference in _PURINES and sample in _PURINES
    both_pyrimidines = reference in _PYRIMIDINES and sample in _PYRIMIDINES
    if both_purines or both_pyrimidines:
        return MutationType.TRANSITION
    return MutationType.TRANSVERSION


def _compare_aligned(reference: str, sample: str) -> list[Variant]:
    """Compare two equal-length (possibly gapped) sequences column by column.

    Args:
        reference: The reference sequence, possibly containing gaps.
        sample: The sample sequence, possibly containing gaps.

    Returns:
        A list of :class:`Variant` for every differing column.
    """
    variants: list[Variant] = []

    for position, (ref_base, sample_base) in enumerate(zip(reference, sample)):
        if ref_base == sample_base:
            continue

        if ref_base == _GAP_CHARACTER:
            mutation_type = MutationType.INSERTION
        elif sample_base == _GAP_CHARACTER:
            mutation_type = MutationType.DELETION
        else:
            mutation_type = _classify_substitution(ref_base, sample_base)

        variants.append(
            Variant(
                position=position,
                reference=ref_base,
                sample=sample_base,
                mutation_type=mutation_type,
            )
        )

    return variants


def compare_sequences(
    reference: str, sample: str, *, force_alignment: bool = False
) -> MutationReport:
    """Compare two DNA sequences and report every difference.

    Sequences of equal length are compared directly. Sequences of differing
    length (or when ``force_alignment`` is set) are aligned first so that
    insertions and deletions are detected rather than shifting every
    downstream base out of register.

    Args:
        reference: The reference DNA sequence.
        sample: The sample DNA sequence to compare against it.
        force_alignment: Align even when the lengths already match.

    Returns:
        A :class:`MutationReport` with all variants and summary statistics.

    Raises:
        InvalidSequenceError: If either sequence is not valid DNA.
    """
    ref = _require_valid_dna(reference)
    smp = _require_valid_dna(sample)

    if force_alignment or len(ref) != len(smp):
        alignment = align(ref, smp)
        ref, smp = alignment.aligned_seq1, alignment.aligned_seq2

    return MutationReport(
        variants=_compare_aligned(ref, smp),
        aligned_reference=ref,
        aligned_sample=smp,
    )


def find_snps(reference: str, sample: str) -> list[Variant]:
    """Return only the single-base substitutions between two sequences.

    Insertions and deletions are excluded; use :func:`compare_sequences` for
    the full picture.

    Args:
        reference: The reference DNA sequence.
        sample: The sample DNA sequence.

    Returns:
        A list of substitution :class:`Variant` only.

    Raises:
        InvalidSequenceError: If either sequence is not valid DNA.
    """
    report = compare_sequences(reference, sample)
    return [
        v
        for v in report.variants
        if v.mutation_type in (MutationType.TRANSITION, MutationType.TRANSVERSION)
    ]
