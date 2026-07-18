"""Unit tests for :mod:`bioseq_toolkit.core.alignment`.

Alignment can have ties (several alignments sharing the optimal score), so
these tests focus on scores and invariant properties rather than demanding
one exact gapped string, except where the result is unambiguous.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.alignment import Alignment, align


class TestIdenticalSequences:
    """Aligning a sequence with itself."""

    @pytest.mark.parametrize("seq", ["ATCG", "A", "ATCGATCGATCG", "GGGG"])
    def test_score_equals_length_for_identical(self, seq: str) -> None:
        # Every column is a match worth +1, so score == length.
        assert align(seq, seq).score == len(seq)

    def test_identical_sequences_have_no_gaps(self) -> None:
        result = align("ATCG", "ATCG")
        assert "-" not in result.aligned_seq1
        assert "-" not in result.aligned_seq2

    def test_identity_is_100_percent(self) -> None:
        assert align("ATCG", "ATCG").identity == pytest.approx(100.0)


class TestSubstitutions:
    """Sequences of equal length differing by substitutions."""

    def test_single_mismatch_score(self) -> None:
        # 3 matches (+3) + 1 mismatch (-1) = 2
        assert align("ATCG", "ATGG").score == 2

    def test_all_mismatches_score(self) -> None:
        # 4 mismatches = -4 (cheaper than gapping everything)
        assert align("AAAA", "CCCC").score == -4

    def test_identity_reflects_mismatch(self) -> None:
        assert align("ATCG", "ATGG").identity == pytest.approx(75.0)


class TestGaps:
    """Sequences of differing length require gaps."""

    def test_deletion_creates_gap(self) -> None:
        result = align("ATCG", "ACG")
        # One sequence must contain exactly one gap character.
        total_gaps = result.aligned_seq1.count("-") + result.aligned_seq2.count("-")
        assert total_gaps == 1

    def test_gap_alignment_score(self) -> None:
        # 3 matches (+3) + 1 gap (-2) = 1
        assert align("ATCG", "ACG").score == 1

    def test_aligned_strings_have_equal_length(self) -> None:
        result = align("ATCGATCG", "ATG")
        assert len(result.aligned_seq1) == len(result.aligned_seq2)


class TestEmptySequences:
    """Edge cases involving empty input."""

    def test_empty_against_sequence_is_all_gaps(self) -> None:
        result = align("", "ATCG")
        assert result.aligned_seq1 == "----"
        assert result.aligned_seq2 == "ATCG"

    def test_empty_against_sequence_score(self) -> None:
        # 4 gaps at -2 each
        assert align("", "ATCG").score == -8

    def test_both_empty(self) -> None:
        result = align("", "")
        assert result.aligned_seq1 == ""
        assert result.score == 0

    def test_empty_identity_is_zero(self) -> None:
        assert align("", "").identity == 0.0


class TestInvariants:
    """Properties that must hold for any alignment."""

    @pytest.mark.parametrize(
        ("s1", "s2"),
        [
            ("ATCG", "ATGG"),
            ("GATTACA", "GCATGCU"),
            ("AAAA", "AA"),
            ("A", "TTTT"),
        ],
    )
    def test_aligned_lengths_always_equal(self, s1: str, s2: str) -> None:
        result = align(s1, s2)
        assert len(result.aligned_seq1) == len(result.aligned_seq2)

    @pytest.mark.parametrize(
        ("s1", "s2"),
        [("ATCG", "ATGG"), ("GATTACA", "GCATGCU"), ("AAAA", "AA")],
    )
    def test_removing_gaps_recovers_originals(self, s1: str, s2: str) -> None:
        """Stripping gap characters must return the input sequences."""
        result = align(s1, s2)
        assert result.aligned_seq1.replace("-", "") == s1
        assert result.aligned_seq2.replace("-", "") == s2

    @pytest.mark.parametrize(
        ("s1", "s2"), [("ATCG", "ATGG"), ("GATTACA", "GCATGCU"), ("AAAA", "AA")]
    )
    def test_score_is_symmetric(self, s1: str, s2: str) -> None:
        """align(a, b) and align(b, a) must produce the same score."""
        assert align(s1, s2).score == align(s2, s1).score


class TestCustomScoring:
    """Tests for user-supplied scoring parameters."""

    def test_custom_match_score(self) -> None:
        assert align("ATCG", "ATCG", match=2).score == 8

    def test_harsher_gap_penalty_still_valid(self) -> None:
        result = align("ATCG", "ACG", gap=-10)
        # With a very harsh gap penalty, mismatches become preferable.
        assert len(result.aligned_seq1) == len(result.aligned_seq2)

    def test_returns_alignment_instance(self) -> None:
        assert isinstance(align("AT", "AT"), Alignment)
