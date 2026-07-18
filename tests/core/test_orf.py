"""Unit tests for :mod:`bioseq_toolkit.core.orf`.

ORF detection has several independent behaviours worth pinning down: correct
coordinates, all three reading frames, both strands, the minimum-length
filter, and correct handling of starts with no downstream stop codon.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.orf import ORF, find_orfs
from bioseq_toolkit.utils.exceptions import InvalidSequenceError

# ATG AAA CGT TAG -> Met-Lys-Arg then stop. 12 bases.
SIMPLE_ORF = "ATGAAACGTTAG"


class TestBasicDetection:
    """Tests for finding a single, simple ORF."""

    def test_finds_simple_orf(self) -> None:
        orfs = find_orfs(SIMPLE_ORF, min_length=6)
        assert len(orfs) == 1

    def test_orf_coordinates(self) -> None:
        orf = find_orfs(SIMPLE_ORF, min_length=6)[0]
        assert orf.start == 0
        assert orf.end == 12
        assert orf.length == 12

    def test_orf_sequence_includes_start_and_stop(self) -> None:
        orf = find_orfs(SIMPLE_ORF, min_length=6)[0]
        assert orf.sequence == SIMPLE_ORF
        assert orf.sequence.startswith("ATG")
        assert orf.sequence.endswith("TAG")

    def test_protein_excludes_stop_codon(self) -> None:
        orf = find_orfs(SIMPLE_ORF, min_length=6)[0]
        assert orf.protein == "MKR"
        assert "*" not in orf.protein

    def test_returns_orf_instances(self) -> None:
        orfs = find_orfs(SIMPLE_ORF, min_length=6)
        assert all(isinstance(o, ORF) for o in orfs)


class TestMinLength:
    """Tests for the minimum-length filter."""

    def test_short_orf_filtered_by_default(self) -> None:
        # The 12-base ORF is below the default 75-base floor.
        assert find_orfs(SIMPLE_ORF) == []

    def test_min_length_boundary_inclusive(self) -> None:
        # min_length == actual length should still be included.
        assert len(find_orfs(SIMPLE_ORF, min_length=12)) == 1

    def test_min_length_just_above_excludes(self) -> None:
        assert find_orfs(SIMPLE_ORF, min_length=13) == []


class TestReadingFrames:
    """Tests that ORFs are found in non-zero reading frames."""

    def test_finds_orf_in_frame_1(self) -> None:
        # One leading base pushes the ORF into frame 1.
        seq = "C" + SIMPLE_ORF
        orfs = find_orfs(seq, min_length=6, both_strands=False)
        assert len(orfs) == 1
        assert orfs[0].frame == 1
        assert orfs[0].start == 1

    def test_finds_orf_in_frame_2(self) -> None:
        seq = "CC" + SIMPLE_ORF
        orfs = find_orfs(seq, min_length=6, both_strands=False)
        assert len(orfs) == 1
        assert orfs[0].frame == 2


class TestStrands:
    """Tests for forward vs reverse strand searching."""

    def test_forward_strand_marked_plus(self) -> None:
        orfs = find_orfs(SIMPLE_ORF, min_length=6, both_strands=False)
        assert orfs[0].strand == "+"

    def test_both_strands_searched_by_default(self) -> None:
        # The reverse complement of this sequence also contains an ORF.
        seq = "CTAACGTTTCAT"  # reverse complement of SIMPLE_ORF
        orfs = find_orfs(seq, min_length=6)
        strands = {o.strand for o in orfs}
        assert "-" in strands

    def test_both_strands_false_searches_forward_only(self) -> None:
        seq = "CTAACGTTTCAT"
        orfs = find_orfs(seq, min_length=6, both_strands=False)
        assert all(o.strand == "+" for o in orfs)


class TestNoOrfCases:
    """Tests for sequences that should yield no ORFs."""

    def test_no_start_codon(self) -> None:
        assert find_orfs("CCCGGGCCCGGG", min_length=3, both_strands=False) == []

    def test_start_without_stop_is_not_an_orf(self) -> None:
        # ATG present but no in-frame stop codon anywhere.
        assert find_orfs("ATGAAACGTAAA", min_length=3, both_strands=False) == []


class TestSorting:
    """Tests that results are ordered longest first."""

    def test_results_sorted_by_length_descending(self) -> None:
        # A long ORF followed by a shorter one.
        seq = "ATGAAACGTACGGGCTTTCCCTAG" + "AAA" + "ATGGGGCCCTAA"
        orfs = find_orfs(seq, min_length=6, both_strands=False)
        lengths = [o.length for o in orfs]
        assert lengths == sorted(lengths, reverse=True)


class TestErrors:
    """Tests for invalid input."""

    @pytest.mark.parametrize("bad", ["", "ATCGX", "AUCG"])
    def test_raises_on_invalid_sequence(self, bad: str) -> None:
        with pytest.raises(InvalidSequenceError):
            find_orfs(bad)
