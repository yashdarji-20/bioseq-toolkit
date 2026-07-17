"""Unit tests for :mod:`bioseq_toolkit.cli.main`.

CLI tests use two pytest tools beyond the usual: ``capsys`` to capture what
the program prints, and ``tmp_path`` to supply a real input file. Because
``main()`` returns an exit code (instead of calling sys.exit directly), we
can assert on that return value directly.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from bioseq_toolkit.cli.main import build_parser, main


@pytest.fixture
def sample_fasta(tmp_path: Path) -> Path:
    """Create a small FASTA file and return its path."""
    f = tmp_path / "sample.fasta"
    f.write_text(">seq1\nATGGCTTAAGGGCCC\n>seq2\nATGAAACGTACGT\n")
    return f


class TestParser:
    """Tests for the argument parser."""

    def test_builds_without_error(self) -> None:
        parser = build_parser()
        assert parser.prog == "bioseq"

    def test_missing_command_exits(self) -> None:
        # argparse exits (SystemExit) when a required subcommand is missing.
        with pytest.raises(SystemExit):
            build_parser().parse_args([])


class TestGcCommand:
    """Tests for the 'gc' subcommand."""

    def test_gc_prints_percentages(
        self, sample_fasta: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        exit_code = main(["gc", str(sample_fasta)])
        out = capsys.readouterr().out
        assert exit_code == 0
        assert "seq1" in out
        assert "GC=60.00%" in out

    def test_gc_returns_zero_on_success(self, sample_fasta: Path) -> None:
        assert main(["gc", str(sample_fasta)]) == 0


class TestTranslateCommand:
    """Tests for the 'translate' subcommand."""

    def test_translate_output(
        self, sample_fasta: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        main(["translate", str(sample_fasta)])
        out = capsys.readouterr().out
        assert "MA" in out  # seq1 translates to MA
        assert "MKRT" in out  # seq2 translates to MKRT


class TestStatsCommand:
    """Tests for the 'stats' subcommand."""

    def test_stats_output(
        self, sample_fasta: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        main(["stats", str(sample_fasta)])
        out = capsys.readouterr().out
        assert "length=15" in out
        assert "A=3" in out


class TestErrorHandling:
    """Tests for graceful error handling."""

    def test_missing_file_returns_error_code(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        exit_code = main(["gc", "does_not_exist.fasta"])
        err = capsys.readouterr().err
        assert exit_code == 1
        assert "not found" in err

    def test_invalid_sequence_returns_error_code(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        bad = tmp_path / "bad.fasta"
        bad.write_text(">seq1\nATCGX\n")  # X is invalid for GC content
        exit_code = main(["gc", str(bad)])
        err = capsys.readouterr().err
        assert exit_code == 1
        assert "Error" in err
