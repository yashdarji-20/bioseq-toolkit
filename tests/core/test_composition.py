"""Unit tests for :mod:`bioseq_toolkit.core.composition`.

Covers GC/AT content and nucleotide frequency, including the two behaviours
that matter most: correct percentages (using ``pytest.approx`` for repeating
decimals) and raising :class:`InvalidSequenceError` on non-DNA input.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.composition import (
    at_content,
    gc_content,
    nucleotide_frequency,
)
from bioseq_toolkit.utils.exceptions import InvalidSequenceError


class TestGcContent:
    """Tests for :func:`gc_content`."""

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("GGCC", 100.0),      # all G/C
            ("AATT", 0.0),        # no G/C
            ("ATCG", 50.0),       # half
            ("  atgc ", 50.0),    # normalized (whitespace + lowercase)
            ("ATGGG", 60.0),      # 3 of 5
            ("ATG", 100 / 3),     # repeating decimal -> needs approx
        ],
    )
    def test_gc_content_value(self, sequence: str, expected: float) -> None:
        assert gc_content(sequence) == pytest.approx(expected)

    @pytest.mark.parametrize("bad", ["", "ATCGX", "AUCG", "1234"])
    def test_gc_content_raises_on_invalid(self, bad: str) -> None:
        with pytest.raises(InvalidSequenceError):
            gc_content(bad)


class TestAtContent:
    """Tests for :func:`at_content`."""

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("AATT", 100.0),
            ("GGCC", 0.0),
            ("ATCG", 50.0),
        ],
    )
    def test_at_content_value(self, sequence: str, expected: float) -> None:
        assert at_content(sequence) == pytest.approx(expected)

    def test_gc_and_at_sum_to_100(self) -> None:
        """For pure A/C/G/T sequences, GC% + AT% should equal 100%."""
        seq = "ATCGATCGGG"
        assert gc_content(seq) + at_content(seq) == pytest.approx(100.0)

    def test_at_content_raises_on_invalid(self) -> None:
        with pytest.raises(InvalidSequenceError):
            at_content("ATCGX")


class TestNucleotideFrequency:
    """Tests for :func:`nucleotide_frequency`."""

    def test_counts_all_bases(self) -> None:
        assert nucleotide_frequency("AATTGC") == {"A": 2, "T": 2, "G": 1, "C": 1}

    def test_absent_bases_are_zero(self) -> None:
        # A, T, C never appear but must still be present with count 0.
        assert nucleotide_frequency("GGG") == {"A": 0, "C": 0, "G": 3, "T": 0}

    def test_always_returns_all_four_keys(self) -> None:
        result = nucleotide_frequency("ACGT")
        assert set(result.keys()) == {"A", "C", "G", "T"}

    def test_raises_on_invalid(self) -> None:
        with pytest.raises(InvalidSequenceError):
            nucleotide_frequency("XYZ")