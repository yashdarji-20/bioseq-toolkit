"""Unit tests for :mod:`bioseq_toolkit.core.validation`.

These tests document the validation contract by example: what counts as a
valid DNA / RNA / protein sequence, and what does not. Edge cases (empty
input, whitespace, case, wrong alphabet) are tested deliberately because
that is where bugs hide.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.validation import (
    is_valid_dna,
    is_valid_protein,
    is_valid_rna,
)


class TestIsValidDna:
    """Tests for :func:`is_valid_dna`."""

    @pytest.mark.parametrize(
        "sequence",
        [
            "ATCG",  # canonical
            "atcg",  # lowercase is normalized
            "  ATCG  ",  # surrounding whitespace is stripped
            "AAAA",  # single-base repeats are fine
            "ATCGATCGATCG",  # longer sequence
        ],
    )
    def test_valid_dna_returns_true(self, sequence: str) -> None:
        assert is_valid_dna(sequence) is True

    @pytest.mark.parametrize(
        "sequence",
        [
            "",  # empty is not a sequence
            "   ",  # whitespace-only is empty after strip
            "ATCGX",  # X is not a base
            "ATCGN",  # N (ambiguity) rejected in v1
            "AUCG",  # U belongs to RNA, not DNA
            "ATCG123",  # digits are not bases
        ],
    )
    def test_invalid_dna_returns_false(self, sequence: str) -> None:
        assert is_valid_dna(sequence) is False


class TestIsValidRna:
    """Tests for :func:`is_valid_rna`."""

    @pytest.mark.parametrize("sequence", ["AUCG", "aucg", "  AUCG\n", "UUUU"])
    def test_valid_rna_returns_true(self, sequence: str) -> None:
        assert is_valid_rna(sequence) is True

    @pytest.mark.parametrize(
        "sequence",
        [
            "",  # empty
            "ATCG",  # T belongs to DNA, not RNA
            "AUCGX",  # X is not a base
        ],
    )
    def test_invalid_rna_returns_false(self, sequence: str) -> None:
        assert is_valid_rna(sequence) is False


class TestIsValidProtein:
    """Tests for :func:`is_valid_protein`."""

    @pytest.mark.parametrize("sequence", ["MKLV", "mklv", "  ACDEFGHIKLMNPQRSTVWY  "])
    def test_valid_protein_returns_true(self, sequence: str) -> None:
        assert is_valid_protein(sequence) is True

    @pytest.mark.parametrize(
        "sequence",
        [
            "",  # empty
            "MKLB",  # B is not one of the 20 standard amino acids
            "MKL1",  # digit
        ],
    )
    def test_invalid_protein_returns_false(self, sequence: str) -> None:
        assert is_valid_protein(sequence) is False
