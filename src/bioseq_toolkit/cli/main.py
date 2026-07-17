"""Command-line interface for BioSeq Toolkit.

This module is intentionally THIN: it parses arguments, reads FASTA input,
delegates to the tested ``core`` functions, formats results, and handles
errors. It contains no sequence-analysis logic of its own.

Usage examples::

    bioseq gc sample.fasta
    bioseq translate sample.fasta
    bioseq reverse sample.fasta
    bioseq stats sample.fasta
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence

from bioseq_toolkit import __version__
from bioseq_toolkit.core.composition import gc_content, nucleotide_frequency
from bioseq_toolkit.core.transforms import reverse_complement
from bioseq_toolkit.core.translation import translate
from bioseq_toolkit.io.fasta import FastaRecord, read_fasta
from bioseq_toolkit.utils.exceptions import BioSeqError
from bioseq_toolkit.utils.logging_config import configure_logging, get_logger

logger = get_logger(__name__)


def _handle_gc(records: list[FastaRecord]) -> None:
    """Print GC content for each record."""
    for record in records:
        print(f"{record.id}\tGC={gc_content(record.sequence):.2f}%")


def _handle_reverse(records: list[FastaRecord]) -> None:
    """Print the reverse complement of each record."""
    for record in records:
        print(f"{record.id}\t{reverse_complement(record.sequence)}")


def _handle_translate(records: list[FastaRecord]) -> None:
    """Print the protein translation of each record."""
    for record in records:
        print(f"{record.id}\t{translate(record.sequence)}")


def _handle_stats(records: list[FastaRecord]) -> None:
    """Print length, GC content, and base counts for each record."""
    for record in records:
        freq = nucleotide_frequency(record.sequence)
        counts = " ".join(f"{base}={freq[base]}" for base in "ACGT")
        print(
            f"{record.id}\tlength={len(record.sequence)}\t"
            f"GC={gc_content(record.sequence):.2f}%\t{counts}"
        )


# Maps each subcommand name to the function that handles it.
_COMMANDS = {
    "gc": _handle_gc,
    "reverse": _handle_reverse,
    "translate": _handle_translate,
    "stats": _handle_stats,
}


def build_parser() -> argparse.ArgumentParser:
    """Build and return the argument parser (separated out for testability)."""
    parser = argparse.ArgumentParser(
        prog="bioseq",
        description="DNA / RNA / protein sequence analysis toolkit.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging."
    )

    subparsers = parser.add_subparsers(dest="command", required=True)
    for name, handler in _COMMANDS.items():
        sub = subparsers.add_parser(name, help=f"Run '{name}' on a FASTA file.")
        sub.add_argument("fasta", help="Path to a FASTA file.")
        sub.set_defaults(handler=handler)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point.

    Args:
        argv: Optional argument list (defaults to ``sys.argv[1:]``). Passing
            it explicitly makes the function testable.

    Returns:
        Process exit code: 0 on success, 1 on a handled BioSeq error.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    configure_logging(logging.DEBUG if args.verbose else logging.INFO)

    try:
        records = read_fasta(args.fasta)
        logger.debug("Read %d record(s) from %s", len(records), args.fasta)
        args.handler(records)
    except FileNotFoundError:
        print(f"Error: file not found: {args.fasta}", file=sys.stderr)
        return 1
    except BioSeqError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
