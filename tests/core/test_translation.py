"""Unit tests for :mod:`bioseq_toolkit.core.translation`.

Each documented design decision gets its own test, so the test file doubles
as executable documentation of how translation is meant to behave.
"""

from __future__ import annotations

import pytest

from bioseq_toolkit.core.translation import translate
from bioseq_toolkit.utils.exceptions import InvalidSequenceError


class TestTranslate:
    """Tests for :func:`translate`."""

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("ATG", "M"),            # single start codon -> Methionine
            ("ATGAAACGT", "MKR"),    # Met-Lys-Arg, a known translation
            ("ATGGCT", "MA"),        # Met-Ala
            ("atgaaacgt", "MKR"),    # lowercase is normalized
        ],
    )
    def test_translate_known_sequences(self, sequence: str, expected: str) -> None:
        assert translate(sequence) == expected

    def test_stops_at_first_stop_codon_by_default(self) -> None:
        # ATG GCT TAA -> M A [stop]; the stop codon ends translation.
        assert translate("ATGGCTTAA") == "MA"

    def test_full_translation_includes_stop_symbol(self) -> None:
        # With to_stop=False the stop codon is represented as "*".
        assert translate("ATGGCTTAA", to_stop=False) == "MA*"

    def test_leading_stop_gives_empty_protein(self) -> None:
        # TAA is a stop codon; with to_stop=True nothing is translated.
        assert translate("TAA") == ""

    def test_leading_stop_full_mode_gives_star(self) -> None:
        assert translate("TAA", to_stop=False) == "*"

    @pytest.mark.parametrize(
        ("sequence", "expected"),
        [
            ("ATGGC", "M"),     # 5 bases: [ATG] translated, "GC" dropped
            ("ATGA", "M"),      # 4 bases: [ATG] translated, "A" dropped
            ("ATGAAAC", "MK"),  # 7 bases: [ATG][AAA] translated, "C" dropped
        ],
    )
    def test_trailing_partial_codon_is_ignored(
        self, sequence: str, expected: str
    ) -> None:
        assert translate(sequence) == expected

    def test_empty_sequence_raises(self) -> None:
        with pytest.raises(InvalidSequenceError):
            translate("")

    def test_invalid_sequence_raises(self) -> None:
        with pytest.raises(InvalidSequenceError):
            translate("ATGX")