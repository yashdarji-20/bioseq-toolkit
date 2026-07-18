"""Unit tests for :mod:`bioseq_toolkit.core.mutation`.

Covers substitution detection, transition/transversion classification, indel
handling via automatic alignment, and the summary statistics. The
classification tests enumerate the full purine/pyrimidine matrix so no
substitution type is left unverified.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.mutation import (
    MutationReport,
    MutationType,
    Variant,
    compare_sequences,
    find_snps,
)
from bioseq_toolkit.utils.exceptions import InvalidSequenceError


class TestIdenticalSequences:
    """Comparing a sequence with itself."""

    def test_no_variants(self) -> None:
        assert compare_sequences("ATCGATCG", "ATCGATCG").total_variants == 0

    def test_identity_is_100(self) -> None:
        assert compare_sequences("ATCG", "ATCG").identity == pytest.approx(100.0)

    def test_variants_list_is_empty(self) -> None:
        assert compare_sequences("ATCG", "ATCG").variants == []


class TestSubstitutionDetection:
    """Finding single-base substitutions."""

    def test_single_substitution_found(self) -> None:
        report = compare_sequences("ATCGATCG", "GTCGATCG")
        assert report.total_variants == 1

    def test_variant_records_position_and_bases(self) -> None:
        variant = compare_sequences("ATCGATCG", "GTCGATCG").variants[0]
        assert variant.position == 0
        assert variant.reference == "A"
        assert variant.sample == "G"

    def test_multiple_substitutions(self) -> None:
        report = compare_sequences("ATCGATCGAT", "GTCGCTCGAT")
        assert [v.position for v in report.variants] == [0, 4]

    def test_variant_str_format(self) -> None:
        variant = compare_sequences("ATCG", "GTCG").variants[0]
        assert str(variant) == "0:A>G"

    def test_identity_reflects_differences(self) -> None:
        # 2 differences in 10 positions -> 80% identity
        report = compare_sequences("ATCGATCGAT", "GTCGCTCGAT")
        assert report.identity == pytest.approx(80.0)


class TestMutationClassification:
    """Transition vs transversion classification."""

    @pytest.mark.parametrize(
        ("reference", "sample"),
        [
            ("A", "G"),  # purine -> purine
            ("G", "A"),
            ("C", "T"),  # pyrimidine -> pyrimidine
            ("T", "C"),
        ],
    )
    def test_transitions(self, reference: str, sample: str) -> None:
        variant = compare_sequences(reference, sample).variants[0]
        assert variant.mutation_type is MutationType.TRANSITION

    @pytest.mark.parametrize(
        ("reference", "sample"),
        [
            ("A", "C"),  # purine -> pyrimidine
            ("A", "T"),
            ("G", "C"),
            ("G", "T"),
            ("C", "A"),  # pyrimidine -> purine
            ("C", "G"),
            ("T", "A"),
            ("T", "G"),
        ],
    )
    def test_transversions(self, reference: str, sample: str) -> None:
        variant = compare_sequences(reference, sample).variants[0]
        assert variant.mutation_type is MutationType.TRANSVERSION


class TestSummaryStatistics:
    """Counts and ratios on the report."""

    def test_counts_transitions_and_transversions(self) -> None:
        report = compare_sequences("ATCGATCGAT", "GTCGCTCGAT")
        assert report.transitions == 1
        assert report.transversions == 1

    def test_ti_tv_ratio(self) -> None:
        report = compare_sequences("ATCGATCGAT", "GTCGCTCGAT")
        assert report.ti_tv_ratio == pytest.approx(1.0)

    def test_ti_tv_ratio_zero_when_no_transversions(self) -> None:
        # Only a transition, so the ratio is defined as 0.0.
        report = compare_sequences("ATCG", "GTCG")
        assert report.ti_tv_ratio == 0.0

    def test_returns_mutation_report(self) -> None:
        assert isinstance(compare_sequences("AT", "AT"), MutationReport)


class TestIndels:
    """Insertions and deletions via automatic alignment."""

    def test_different_lengths_are_aligned(self) -> None:
        report = compare_sequences("ATCGATCG", "ATCGTCG")
        # Alignment pads the shorter sequence, so both are equal length.
        assert len(report.aligned_reference) == len(report.aligned_sample)

    def test_deletion_detected(self) -> None:
        report = compare_sequences("ATCGATCG", "ATCGTCG")
        assert report.indels == 1
        assert report.variants[0].mutation_type is MutationType.DELETION

    def test_insertion_detected(self) -> None:
        # Sample is longer, so the reference gains a gap.
        report = compare_sequences("ATCGTCG", "ATCGATCG")
        assert report.indels == 1
        assert report.variants[0].mutation_type is MutationType.INSERTION

    def test_deletion_not_reported_as_many_substitutions(self) -> None:
        """A single deletion must not cascade into false SNPs."""
        report = compare_sequences("ATCGATCG", "ATCGTCG")
        assert report.total_variants == 1

    def test_force_alignment_on_equal_lengths(self) -> None:
        report = compare_sequences("ATCG", "ATCG", force_alignment=True)
        assert report.total_variants == 0


class TestFindSnps:
    """The substitution-only helper."""

    def test_returns_only_substitutions(self) -> None:
        snps = find_snps("ATCGATCGAT", "GTCGCTCGAT")
        assert len(snps) == 2

    def test_excludes_indels(self) -> None:
        # The only difference here is a deletion, so no SNPs.
        assert find_snps("ATCGATCG", "ATCGTCG") == []

    def test_returns_variant_instances(self) -> None:
        assert all(isinstance(s, Variant) for s in find_snps("ATCG", "GTCG"))


class TestErrors:
    """Invalid input handling."""

    @pytest.mark.parametrize("bad", ["", "ATCGX", "AUCG"])
    def test_invalid_reference_raises(self, bad: str) -> None:
        with pytest.raises(InvalidSequenceError):
            compare_sequences(bad, "ATCG")

    @pytest.mark.parametrize("bad", ["", "ATCGX", "AUCG"])
    def test_invalid_sample_raises(self, bad: str) -> None:
        with pytest.raises(InvalidSequenceError):
            compare_sequences("ATCG", bad)
