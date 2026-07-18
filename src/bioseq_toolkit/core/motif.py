"""Motif and restriction-site searching for BioSeq Toolkit.

Both features reduce to the same problem: find every position where a pattern
occurs in a sequence, where the pattern may contain IUPAC ambiguity codes
(``N`` = any base, ``R`` = A/G, and so on). :func:`find_pattern` implements
that once; motif finding and restriction-enzyme searching both build on it.

Overlapping occurrences are reported. Searching ``"AAAA"`` for ``"AA"`` yields
matches at positions 0, 1 and 2, because overlapping sites are biologically
real and a naive non-overlapping scan would miss them.
"""

from __future__ import annotations

from dataclasses import dataclass

from bioseq_toolkit.core.composition import _require_valid_dna

# IUPAC nucleotide ambiguity codes: each maps to the set of bases it matches.
IUPAC_CODES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),  # puRine
    "Y": frozenset("CT"),  # pYrimidine
    "S": frozenset("GC"),  # Strong (3 H-bonds)
    "W": frozenset("AT"),  # Weak (2 H-bonds)
    "K": frozenset("GT"),  # Keto
    "M": frozenset("AC"),  # aMino
    "B": frozenset("CGT"),  # not A
    "D": frozenset("AGT"),  # not C
    "H": frozenset("ACT"),  # not G
    "V": frozenset("ACG"),  # not T
    "N": frozenset("ACGT"),  # aNy base
}

# A small database of common restriction enzymes.
# Each entry maps the enzyme name to (recognition site, cut offset), where the
# offset is how many bases from the start of the site the cut occurs.
RESTRICTION_ENZYMES: dict[str, tuple[str, int]] = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "NotI": ("GCGGCCGC", 2),
    "PstI": ("CTGCAG", 5),
    "SmaI": ("CCCGGG", 3),
    "EcoRV": ("GATATC", 3),
    "HaeIII": ("GGCC", 2),
}


@dataclass(frozen=True)
class Match:
    """A single pattern occurrence.

    Attributes:
        start: 0-based index where the match begins.
        end: 0-based exclusive index where the match ends.
        matched: The actual sequence text that matched.
    """

    start: int
    end: int
    matched: str


def _matches_at(sequence: str, pattern: str, index: int) -> bool:
    """Return ``True`` if ``pattern`` matches ``sequence`` starting at ``index``."""
    for offset, code in enumerate(pattern):
        allowed = IUPAC_CODES.get(code)
        if allowed is None or sequence[index + offset] not in allowed:
            return False
    return True


def find_pattern(sequence: str, pattern: str) -> list[Match]:
    """Find every (possibly overlapping) occurrence of a pattern.

    The pattern may contain IUPAC ambiguity codes such as ``N`` or ``R``.

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).
        pattern: The pattern to search for, in IUPAC codes.

    Returns:
        A list of :class:`Match`, in ascending position order.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
        ValueError: If ``pattern`` is empty or contains an unknown code.
    """
    seq = _require_valid_dna(sequence)
    pat = pattern.strip().upper()

    if not pat:
        raise ValueError("Pattern must not be empty.")
    unknown = set(pat) - set(IUPAC_CODES)
    if unknown:
        raise ValueError(f"Unknown IUPAC code(s) in pattern: {sorted(unknown)}")

    matches: list[Match] = []
    # Step by 1 (not by len(pat)) so overlapping occurrences are all found.
    for i in range(len(seq) - len(pat) + 1):
        if _matches_at(seq, pat, i):
            matches.append(Match(i, i + len(pat), seq[i : i + len(pat)]))
    return matches


def find_motif(sequence: str, motif: str) -> list[Match]:
    """Find all occurrences of a motif in a sequence.

    A thin, well-named alias for :func:`find_pattern`, kept separate because
    "motif" is the term biologists use for this operation.

    Args:
        sequence: A DNA sequence.
        motif: The motif to search for (IUPAC codes allowed).

    Returns:
        A list of :class:`Match`.
    """
    return find_pattern(sequence, motif)


def find_restriction_sites(sequence: str, enzyme: str) -> list[Match]:
    """Find all recognition sites for a named restriction enzyme.

    Args:
        sequence: A DNA sequence.
        enzyme: Enzyme name, e.g. ``"EcoRI"`` (see :data:`RESTRICTION_ENZYMES`).

    Returns:
        A list of :class:`Match` marking each recognition site.

    Raises:
        ValueError: If the enzyme name is not in the database.
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    if enzyme not in RESTRICTION_ENZYMES:
        raise ValueError(
            f"Unknown enzyme {enzyme!r}. " f"Available: {sorted(RESTRICTION_ENZYMES)}"
        )
    site, _cut_offset = RESTRICTION_ENZYMES[enzyme]
    return find_pattern(sequence, site)


def digest(sequence: str, enzyme: str) -> list[str]:
    """Return the DNA fragments produced by cutting with a restriction enzyme.

    Args:
        sequence: A DNA sequence.
        enzyme: Enzyme name, e.g. ``"EcoRI"``.

    Returns:
        The fragments in order. A sequence with no recognition site yields a
        single fragment: the whole sequence.

    Raises:
        ValueError: If the enzyme name is not in the database.
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)
    sites = find_restriction_sites(seq, enzyme)
    _site, cut_offset = RESTRICTION_ENZYMES[enzyme]

    fragments: list[str] = []
    previous_cut = 0
    for match in sites:
        cut_position = match.start + cut_offset
        fragments.append(seq[previous_cut:cut_position])
        previous_cut = cut_position
    fragments.append(seq[previous_cut:])
    return fragments
