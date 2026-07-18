"""Unit tests for :mod:`bioseq_toolkit.core.primer`.

Includes a regression test for a real bug: an earlier implementation of
self-complementarity compared every region against every other, so a region
matched itself trivially and no primer could ever pass. The "known-good
primer must be suitable" test below is what catches that class of error.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.primer import (
    MAX_HOMOPOLYMER_RUN,
    PrimerAnalysis,
    _longest_homopolymer_run,
    _self_complementarity,
    analyze_primer,
)
from bioseq_toolkit.utils.exceptions import InvalidSequenceError

# A primer meeting every design criterion: 20 bases, 50% GC, Tm 60 C,
# 3' GC clamp, no long runs, no self-complementarity.
GOOD_PRIMER = "ACGTTGCATTGACCTGAAGC"


class TestGoodPrimer:
    """A primer that should pass every check."""

    def test_known_good_primer_is_suitable(self) -> None:
        """Regression test: a valid primer must not be rejected."""
        assert analyze_primer(GOOD_PRIMER).is_suitable is True

    def test_good_primer_has_no_warnings(self) -> None:
        assert analyze_primer(GOOD_PRIMER).warnings == []

    def test_measurements_are_recorded(self) -> None:
        result = analyze_primer(GOOD_PRIMER)
        assert result.length == 20
        assert result.gc_percent == pytest.approx(50.0)
        assert result.has_gc_clamp is True

    def test_returns_primer_analysis(self) -> None:
        assert isinstance(analyze_primer(GOOD_PRIMER), PrimerAnalysis)


class TestLengthChecks:
    """Length must fall in the 18-25 base window."""

    def test_too_short_is_flagged(self) -> None:
        result = analyze_primer("ATCGATCGGC")  # 10 bases
        assert any("Too short" in w for w in result.warnings)

    def test_too_long_is_flagged(self) -> None:
        result = analyze_primer("ACGTTGCATTGACCTGAAGCACGTTGCATTGC")  # 32 bases
        assert any("Too long" in w for w in result.warnings)

    def test_length_recorded_correctly(self) -> None:
        assert analyze_primer("ATCGATCGGC").length == 10


class TestGcChecks:
    """GC content must fall in the 40-60% window."""

    def test_low_gc_is_flagged(self) -> None:
        result = analyze_primer("ATATATATATATATATATAT")
        assert any("GC content too low" in w for w in result.warnings)

    def test_high_gc_is_flagged(self) -> None:
        result = analyze_primer("GCGCGCGCGCGCGCGCGCGC")
        assert any("GC content too high" in w for w in result.warnings)


class TestGcClamp:
    """The 3' end should be G or C."""

    def test_missing_clamp_is_flagged(self) -> None:
        result = analyze_primer("ACGTTGCATTGACCTGAAGA")  # ends in A
        assert result.has_gc_clamp is False
        assert any("GC clamp" in w for w in result.warnings)

    @pytest.mark.parametrize("last_base", ["G", "C"])
    def test_clamp_present(self, last_base: str) -> None:
        result = analyze_primer("ACGTTGCATTGACCTGAAG" + last_base)
        assert result.has_gc_clamp is True


class TestHomopolymerRuns:
    """Long single-base runs cause mispriming."""

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("ATCG", 1),
            ("AATCG", 2),
            ("AAACGT", 3),
            ("ATCGAAAAAAA", 7),
            ("GGGG", 4),
        ],
    )
    def test_longest_run_calculation(self, sequence: str, expected: int) -> None:
        assert _longest_homopolymer_run(sequence) == expected

    def test_empty_sequence_run_is_zero(self) -> None:
        assert _longest_homopolymer_run("") == 0

    def test_long_run_is_flagged(self) -> None:
        result = analyze_primer("ATCGAAAAAAATCGATCGGC")
        assert result.max_homopolymer_run > MAX_HOMOPOLYMER_RUN
        assert any("Homopolymer run" in w for w in result.warnings)


class TestSelfComplementarity:
    """Hairpin and primer-dimer detection."""

    def test_good_primer_has_low_self_complementarity(self) -> None:
        """The bug fix: a normal primer must score low, not maximal."""
        assert _self_complementarity(GOOD_PRIMER) <= 4

    def test_hairpin_is_detected(self) -> None:
        # 3' end ATCGATCG is complementary to the 5' start, with a TTTT loop.
        assert _self_complementarity("ATCGATCGTTTTCGATCGAT") > 4

    def test_hairpin_primer_is_flagged(self) -> None:
        result = analyze_primer("ATCGATCGTTTTCGATCGAT")
        assert any("Self-complementary" in w for w in result.warnings)

    def test_single_base_has_no_self_complementarity(self) -> None:
        assert _self_complementarity("A") == 0


class TestReport:
    """The human-readable summary."""

    def test_report_shows_suitable_verdict(self) -> None:
        assert "SUITABLE" in analyze_primer(GOOD_PRIMER).report()

    def test_report_shows_rejected_verdict(self) -> None:
        assert "NOT RECOMMENDED" in analyze_primer("ATATATATATATATATATAT").report()

    def test_report_lists_warnings(self) -> None:
        report = analyze_primer("ATATATATATATATATATAT").report()
        assert "Warnings:" in report


class TestErrors:
    """Invalid input still raises."""

    @pytest.mark.parametrize("bad", ["", "ATCGX", "AUCG", "1234"])
    def test_invalid_sequence_raises(self, bad: str) -> None:
        with pytest.raises(InvalidSequenceError):
            analyze_primer(bad)
