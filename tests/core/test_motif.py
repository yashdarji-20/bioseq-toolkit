"""Unit tests for :mod:`bioseq_toolkit.core.motif`.

Covers pattern matching with IUPAC ambiguity codes, overlapping occurrences,
restriction-site detection, and digestion. The digest tests include the key
invariant: rejoining all fragments must reproduce the original sequence.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.motif import (
    RESTRICTION_ENZYMES,
    Match,
    digest,
    find_motif,
    find_pattern,
    find_restriction_sites,
)
from bioseq_toolkit.utils.exceptions import InvalidSequenceError


class TestExactMatching:
    """Pattern matching with no ambiguity codes."""

    def test_finds_single_occurrence(self) -> None:
        matches = find_pattern("AAATCGAAA", "ATCG")
        assert len(matches) == 1
        assert matches[0].start == 2

    def test_finds_multiple_occurrences(self) -> None:
        matches = find_pattern("ATCGATCGATCG", "ATCG")
        assert [m.start for m in matches] == [0, 4, 8]

    def test_no_match_returns_empty_list(self) -> None:
        assert find_pattern("AAAACCCC", "GGGG") == []

    def test_match_records_correct_span(self) -> None:
        match = find_pattern("AAATCG", "ATCG")[0]
        assert match.start == 2
        assert match.end == 6
        assert match.matched == "ATCG"

    def test_returns_match_instances(self) -> None:
        assert all(isinstance(m, Match) for m in find_pattern("ATCG", "AT"))

    def test_pattern_longer_than_sequence(self) -> None:
        assert find_pattern("AT", "ATCGATCG") == []

    def test_case_insensitive(self) -> None:
        assert len(find_pattern("atcgatcg", "atcg")) == 2


class TestOverlappingMatches:
    """Overlapping occurrences must all be reported."""

    def test_overlapping_matches_found(self) -> None:
        # "AAAA" contains "AA" at positions 0, 1 and 2.
        assert [m.start for m in find_pattern("AAAA", "AA")] == [0, 1, 2]

    def test_overlapping_longer_pattern(self) -> None:
        assert [m.start for m in find_pattern("AAAAA", "AAA")] == [0, 1, 2]


class TestIupacCodes:
    """Pattern matching with ambiguity codes."""

    def test_n_matches_any_base(self) -> None:
        matches = find_pattern("GGATTAGGCTTA", "GGNTTA")
        assert [m.matched for m in matches] == ["GGATTA", "GGCTTA"]

    def test_r_matches_purines_only(self) -> None:
        # R = A or G, so "AR" matches "AA" and "AG" but not "AT" or "AC".
        matches = find_pattern("ATAGACAA", "AR")
        assert all(m.matched[1] in "AG" for m in matches)

    def test_y_matches_pyrimidines_only(self) -> None:
        matches = find_pattern("ACATAG", "AY")
        assert all(m.matched[1] in "CT" for m in matches)

    def test_n_only_pattern_matches_everywhere(self) -> None:
        # "N" matches every single base.
        assert len(find_pattern("ATCG", "N")) == 4


class TestFindMotifAlias:
    """find_motif should behave identically to find_pattern."""

    def test_alias_matches_find_pattern(self) -> None:
        assert find_motif("ATCGATCG", "ATCG") == find_pattern("ATCGATCG", "ATCG")


class TestPatternErrors:
    """Invalid patterns and sequences."""

    def test_empty_pattern_raises(self) -> None:
        with pytest.raises(ValueError):
            find_pattern("ATCG", "")

    def test_unknown_iupac_code_raises(self) -> None:
        with pytest.raises(ValueError):
            find_pattern("ATCG", "ATZQ")

    def test_invalid_sequence_raises(self) -> None:
        with pytest.raises(InvalidSequenceError):
            find_pattern("ATCGX", "AT")


class TestRestrictionSites:
    """Finding recognition sites for named enzymes."""

    def test_finds_ecori_site(self) -> None:
        matches = find_restriction_sites("AAAGAATTCTTT", "EcoRI")
        assert len(matches) == 1
        assert matches[0].matched == "GAATTC"
        assert matches[0].start == 3

    def test_finds_multiple_sites(self) -> None:
        seq = "AAAGAATTCGGGGAATTCTTT"
        assert len(find_restriction_sites(seq, "EcoRI")) == 2

    def test_no_sites_returns_empty(self) -> None:
        assert find_restriction_sites("AAACCCGGG", "EcoRI") == []

    def test_unknown_enzyme_raises(self) -> None:
        with pytest.raises(ValueError):
            find_restriction_sites("ATCG", "NotARealEnzyme")

    @pytest.mark.parametrize("enzyme", list(RESTRICTION_ENZYMES))
    def test_every_enzyme_finds_its_own_site(self, enzyme: str) -> None:
        """Each enzyme must recognise a sequence containing its own site."""
        site, _offset = RESTRICTION_ENZYMES[enzyme]
        matches = find_restriction_sites("AAA" + site + "TTT", enzyme)
        assert len(matches) == 1


class TestDigest:
    """Cutting a sequence into fragments."""

    def test_single_cut_gives_two_fragments(self) -> None:
        assert len(digest("AAAGAATTCTTT", "EcoRI")) == 2

    def test_two_cuts_give_three_fragments(self) -> None:
        assert len(digest("AAAGAATTCGGGGAATTCTTT", "EcoRI")) == 3

    def test_no_site_gives_whole_sequence(self) -> None:
        assert digest("AAACCCGGG", "EcoRI") == ["AAACCCGGG"]

    def test_ecori_cuts_at_correct_offset(self) -> None:
        # EcoRI recognises GAATTC and cuts after the G (offset 1).
        assert digest("AAAGAATTCGGG", "EcoRI") == ["AAAG", "AATTCGGG"]

    @pytest.mark.parametrize(
        "sequence",
        [
            "AAAGAATTCTTT",
            "AAAGAATTCGGGGAATTCTTT",
            "GAATTCGAATTCGAATTC",
            "AAACCCGGG",
        ],
    )
    def test_fragments_rejoin_to_original(self, sequence: str) -> None:
        """The key invariant: concatenating fragments recovers the input."""
        assert "".join(digest(sequence, "EcoRI")) == sequence

    def test_unknown_enzyme_raises(self) -> None:
        with pytest.raises(ValueError):
            digest("ATCG", "NotARealEnzyme")
