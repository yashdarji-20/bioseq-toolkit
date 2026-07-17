"""Unit tests for :mod:`bioseq_toolkit.core.stats`.

Covers the Wallace-rule melting temperature, molecular-weight estimation
(including the single-base edge case), and the SequenceStats bundle.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.stats import (
    SequenceStats,
    melting_temperature,
    molecular_weight,
    sequence_stats,
)
from bioseq_toolkit.utils.exceptions import InvalidSequenceError


class TestMeltingTemperature:
    """Tests for :func:`melting_temperature` (Wallace rule)."""

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("ATGC", 12),  # 2*(A+T=2) + 4*(G+C=2) = 4 + 8
            ("ATATAT", 12),  # 2*6, no G/C
            ("GGGGCC", 24),  # 4*6, no A/T
            ("atgc", 12),  # normalized
        ],
    )
    def test_tm_value(self, sequence: str, expected: int) -> None:
        assert melting_temperature(sequence) == expected

    def test_tm_raises_on_invalid(self) -> None:
        with pytest.raises(InvalidSequenceError):
            melting_temperature("ATCGX")


class TestMolecularWeight:
    """Tests for :func:`molecular_weight`."""

    def test_single_base_has_no_bond_subtraction(self) -> None:
        # A 1-base sequence forms 0 bonds, so no water is subtracted.
        assert molecular_weight("A") == pytest.approx(331.2222)

    def test_weight_increases_with_length(self) -> None:
        assert molecular_weight("ATGC") > molecular_weight("ATG")

    def test_weight_is_positive(self) -> None:
        assert molecular_weight("ATGCATGCATGC") > 0

    def test_raises_on_invalid(self) -> None:
        with pytest.raises(InvalidSequenceError):
            molecular_weight("")


class TestSequenceStats:
    """Tests for :func:`sequence_stats`."""

    def test_returns_sequence_stats_instance(self) -> None:
        result = sequence_stats("ATGC")
        assert isinstance(result, SequenceStats)

    def test_stats_values(self) -> None:
        result = sequence_stats("ATGGCTTAAGGGCCC")
        assert result.length == 15
        assert result.gc_percent == pytest.approx(60.0)
        assert result.at_percent == pytest.approx(40.0)
        assert result.base_counts == {"A": 3, "C": 4, "G": 5, "T": 3}

    def test_gc_at_sum_to_100(self) -> None:
        result = sequence_stats("ATCGATCGGG")
        assert result.gc_percent + result.at_percent == pytest.approx(100.0)

    def test_stats_is_frozen(self) -> None:
        result = sequence_stats("ATGC")
        with pytest.raises(Exception):
            result.length = 99  # type: ignore[misc]

    def test_raises_on_invalid(self) -> None:
        with pytest.raises(InvalidSequenceError):
            sequence_stats("XYZ")
