"""FASTA file reading and writing for BioSeq Toolkit.

FASTA is the standard plain-text format for biological sequences::

    >identifier optional description
    ATCGATCG
    ATCG

This module parses FASTA into :class:`FastaRecord` objects and writes them
back out. It deliberately does NOT validate the biology of sequences: parsing
and validation are separate concerns, and a FASTA file may hold DNA, RNA, or
protein. Validate later with the ``core`` functions if needed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bioseq_toolkit.utils.exceptions import FileFormatError

# How many characters per line when writing sequences (a common FASTA default).
_DEFAULT_LINE_WIDTH = 60


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record: an identifier and its sequence.

    Attributes:
        id: The record identifier (the text after ``>``, first word onward).
        sequence: The sequence with all line breaks removed.
    """

    id: str
    sequence: str


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read a FASTA file into a list of records.

    Multi-line sequences are joined; blank lines are ignored.

    Args:
        path: Path to a FASTA file.

    Returns:
        A list of :class:`FastaRecord`, in file order.

    Raises:
        FileFormatError: If sequence data appears before any header, or a
            header has no sequence, or the file contains no records.
    """
    records: list[FastaRecord] = []
    current_id: str | None = None
    current_seq: list[str] = []

    def _flush() -> None:
        """Finalize the record being accumulated, if any."""
        if current_id is not None:
            if not current_seq:
                raise FileFormatError(f"Header '{current_id}' has no sequence.")
            records.append(FastaRecord(current_id, "".join(current_seq)))

    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue  # skip blank lines
            if line.startswith(">"):
                _flush()                       # close the previous record
                current_id = line[1:].strip()  # drop '>' and surrounding space
                current_seq = []
            else:
                if current_id is None:
                    raise FileFormatError("Sequence data found before any header.")
                current_seq.append(line)

    _flush()  # close the final record

    if not records:
        raise FileFormatError("No FASTA records found in file.")
    return records


def write_fasta(
    records: list[FastaRecord],
    path: str | Path,
    *,
    line_width: int = _DEFAULT_LINE_WIDTH,
) -> None:
    """Write records to a FASTA file, wrapping sequences to a fixed width.

    Args:
        records: The records to write.
        path: Destination file path.
        line_width: Maximum characters per sequence line (default 60).
    """
    with open(path, "w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.id}\n")
            seq = record.sequence
            for i in range(0, len(seq), line_width):
                handle.write(seq[i : i + line_width] + "\n")