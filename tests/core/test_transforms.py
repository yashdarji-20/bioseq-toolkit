"""Unit tests for :mod:`bioseq_toolkit.core.transforms`.

Beyond checking known input/output pairs, these tests verify a mathematical
invariant: applying the reverse complement twice returns the original
sequence. Invariant tests catch whole classes of bugs that specific-value
tests can miss.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.transforms import (
    complement,
    reverse_complement,
    transcribe,
)
from bioseq_toolkit.utils.exceptions import InvalidSequenceError


class TestComplement:
    """Tests for :func:`complement`."""

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("ATCG", "TAGC"),
            ("AAAA", "TTTT"),
            ("GGCC", "CCGG"),
            ("atcg", "TAGC"),  # lowercase normalized
        ],
    )
    def test_complement_value(self, sequence: str, expected: str) -> None:
        assert complement(sequence) == expected

    def test_complement_raises_on_invalid(self) -> None:
        with pytest.raises(InvalidSequenceError):
            complement("ATCGX")


class TestReverseComplement:
    """Tests for :func:`reverse_complement`."""

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("ATCG", "CGAT"),
            ("AAAACCCGGT", "ACCGGGTTTT"),
            ("GGGG", "CCCC"),
            ("  atcg ", "CGAT"),  # normalized
        ],
    )
    def test_reverse_complement_value(self, sequence: str, expected: str) -> None:
        assert reverse_complement(sequence) == expected

    @pytest.mark.parametrize(
        "sequence",
        ["ATCG", "AAAACCCGGT", "GCGCGCTATA", "A", "ACGTACGTACGT"],
    )
    def test_double_reverse_complement_is_identity(self, sequence: str) -> None:
        """revcomp(revcomp(x)) must return the original x."""
        assert reverse_complement(reverse_complement(sequence)) == sequence

    def test_reverse_complement_raises_on_invalid(self) -> None:
        with pytest.raises(InvalidSequenceError):
            reverse_complement("")


class TestTranscribe:
    """Tests for :func:`transcribe`."""

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("ATCG", "AUCG"),
            ("TTTT", "UUUU"),
            ("ACGACG", "ACGACG"),  # no T -> unchanged
            ("atcg", "AUCG"),  # normalized
        ],
    )
    def test_transcribe_value(self, sequence: str, expected: str) -> None:
        assert transcribe(sequence) == expected

    def test_transcribe_has_no_t(self) -> None:
        """Transcribed RNA must never contain the DNA base T."""
        assert "T" not in transcribe("ATCGATCGTTTT")

    def test_transcribe_raises_on_invalid(self) -> None:
        with pytest.raises(InvalidSequenceError):
            transcribe("XYZ")
