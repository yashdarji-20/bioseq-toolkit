"""Unit tests for :mod:`bioseq_toolkit.io.fasta`.

Uses pytest's ``tmp_path`` fixture so each test gets a private temporary
directory. No test files are committed to the repo, and nothing is left
behind after the tests run.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from bioseq_toolkit.io.fasta import FastaRecord, read_fasta, write_fasta
from bioseq_toolkit.utils.exceptions import BioSeqError, FileFormatError


class TestReadFasta:
    """Tests for :func:`read_fasta`."""

    def test_reads_single_record(self, tmp_path: Path) -> None:
        f = tmp_path / "one.fasta"
        f.write_text(">seq1 desc\nATCGATCG\n")
        records = read_fasta(f)
        assert records == [FastaRecord("seq1 desc", "ATCGATCG")]

    def test_reads_multiple_records(self, tmp_path: Path) -> None:
        f = tmp_path / "many.fasta"
        f.write_text(">a\nAAAA\n>b\nCCCC\n")
        records = read_fasta(f)
        assert [r.id for r in records] == ["a", "b"]
        assert [r.sequence for r in records] == ["AAAA", "CCCC"]

    def test_joins_multiline_sequence(self, tmp_path: Path) -> None:
        f = tmp_path / "multiline.fasta"
        f.write_text(">seq1\nATCG\nATCG\nGGGG\n")
        records = read_fasta(f)
        assert records[0].sequence == "ATCGATCGGGGG"

    def test_ignores_blank_lines(self, tmp_path: Path) -> None:
        f = tmp_path / "blanks.fasta"
        f.write_text(">seq1\n\nATCG\n\n\nGGGG\n")
        records = read_fasta(f)
        assert records[0].sequence == "ATCGGGGG"

    def test_sequence_before_header_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "bad1.fasta"
        f.write_text("ATCG\n>seq1\nGGGG\n")
        with pytest.raises(FileFormatError):
            read_fasta(f)

    def test_header_without_sequence_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "bad2.fasta"
        f.write_text(">seq1\n")
        with pytest.raises(FileFormatError):
            read_fasta(f)

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        f = tmp_path / "empty.fasta"
        f.write_text("")
        with pytest.raises(FileFormatError):
            read_fasta(f)

    def test_fileformaterror_is_bioseqerror(self, tmp_path: Path) -> None:
        """The specific error is catchable via the library's base class."""
        f = tmp_path / "empty.fasta"
        f.write_text("")
        with pytest.raises(BioSeqError):
            read_fasta(f)


class TestWriteFasta:
    """Tests for :func:`write_fasta`."""

    def test_write_then_read_round_trip(self, tmp_path: Path) -> None:
        records = [FastaRecord("s1", "ATCGATCG"), FastaRecord("s2", "GGGGCCCC")]
        f = tmp_path / "out.fasta"
        write_fasta(records, f)
        assert read_fasta(f) == records

    def test_wraps_long_sequences(self, tmp_path: Path) -> None:
        # A 10-char sequence with line_width=4 should wrap into 3 lines.
        records = [FastaRecord("s1", "ATCGATCGAT")]
        f = tmp_path / "wrapped.fasta"
        write_fasta(records, f, line_width=4)
        lines = f.read_text().splitlines()
        assert lines == [">s1", "ATCG", "ATCG", "AT"]

    def test_round_trip_preserves_sequence_despite_wrapping(
        self, tmp_path: Path
    ) -> None:
        """Wrapping changes line layout but not the joined sequence."""
        original = FastaRecord("s1", "A" * 150)
        f = tmp_path / "long.fasta"
        write_fasta([original], f, line_width=60)
        assert read_fasta(f)[0].sequence == "A" * 150