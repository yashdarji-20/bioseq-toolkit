"""PCR primer analysis for BioSeq Toolkit.

Evaluates whether a DNA sequence is suitable as a PCR primer by checking the
properties that determine whether amplification will actually work: length,
GC content, melting temperature, the presence of a 3' GC clamp, self-
complementarity (hairpins and primer-dimers), and homopolymer runs.

This module is mostly *orchestration*: it reuses the tested composition,
stats, and transform functions rather than reimplementing any of them.

A suboptimal primer is not an error, so problems are reported as warnings on
the returned :class:`PrimerAnalysis` rather than raised as exceptions. Only
genuinely invalid DNA raises.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from bioseq_toolkit.core.composition import _require_valid_dna, gc_content
from bioseq_toolkit.core.stats import melting_temperature
from bioseq_toolkit.core.transforms import reverse_complement

# Accepted ranges for a usable PCR primer. Collected here so requirements are
# visible and adjustable in one place rather than buried in conditionals.
MIN_LENGTH = 18
MAX_LENGTH = 25
MIN_GC_PERCENT = 40.0
MAX_GC_PERCENT = 60.0
MIN_MELTING_TEMP = 50.0
MAX_MELTING_TEMP = 65.0

# Longest run of one repeated base that is still acceptable.
MAX_HOMOPOLYMER_RUN = 4

# Longest self-complementary stretch that is still acceptable.
MAX_SELF_COMPLEMENTARITY = 4

_GC_BASES = frozenset("GC")


@dataclass(frozen=True)
class PrimerAnalysis:
    """The result of evaluating a candidate PCR primer.

    Attributes:
        sequence: The normalized primer sequence.
        length: Primer length in bases.
        gc_percent: GC content as a percentage.
        melting_temp: Estimated melting temperature in degrees Celsius.
        has_gc_clamp: Whether the 3' (last) base is G or C.
        max_homopolymer_run: Length of the longest single-base run.
        self_complementarity: Length of the longest self-complementary stretch.
        warnings: Human-readable descriptions of every problem found.
    """

    sequence: str
    length: int
    gc_percent: float
    melting_temp: float
    has_gc_clamp: bool
    max_homopolymer_run: int
    self_complementarity: int
    warnings: list[str] = field(default_factory=list)

    @property
    def is_suitable(self) -> bool:
        """Whether the primer passed every check."""
        return not self.warnings

    def report(self) -> str:
        """Return a multi-line human-readable summary."""
        verdict = "SUITABLE" if self.is_suitable else "NOT RECOMMENDED"
        lines = [
            f"Primer: {self.sequence}",
            f"Verdict: {verdict}",
            f"  Length:            {self.length} bases",
            f"  GC content:        {self.gc_percent:.1f}%",
            f"  Melting temp:      {self.melting_temp:.1f} C",
            f"  GC clamp (3' end): {'yes' if self.has_gc_clamp else 'no'}",
            f"  Longest base run:  {self.max_homopolymer_run}",
            f"  Self-complement:   {self.self_complementarity}",
        ]
        if self.warnings:
            lines.append("  Warnings:")
            lines.extend(f"    - {w}" for w in self.warnings)
        return "\n".join(lines)


def _longest_homopolymer_run(sequence: str) -> int:
    """Return the length of the longest run of a single repeated base.

    Args:
        sequence: A normalized DNA sequence.

    Returns:
        The longest run length, or 0 for an empty sequence.
    """
    if not sequence:
        return 0

    longest = 1
    current = 1
    for previous, base in zip(sequence, sequence[1:]):
        current = current + 1 if base == previous else 1
        longest = max(longest, current)
    return longest


def _self_complementarity(sequence: str) -> int:
    """Return the length of the longest 3'-end self-complementary stretch.

    A primer whose 3' end is complementary to an earlier part of itself can
    fold into a hairpin or bind a second copy of itself (a primer-dimer),
    either of which stops it binding the intended target.

    Only the 3' end is checked, and only against regions that do not overlap
    it, because a region is trivially complementary to itself.

    Args:
        sequence: A normalized DNA sequence.

    Returns:
        The length of the longest 3'-end stretch whose complement appears in
        the non-overlapping upstream portion, or 0 if none does.
    """
    length = len(sequence)
    # Longest first: report the worst case found.
    for size in range(length // 2, 0, -1):
        three_prime_end = sequence[length - size :]
        # The complement this end would pair with.
        partner = reverse_complement(three_prime_end)
        # Search only the upstream region that does not overlap the 3' end.
        upstream = sequence[: length - size]
        if partner in upstream:
            return size
    return 0


def analyze_primer(sequence: str) -> PrimerAnalysis:
    """Evaluate a candidate PCR primer against standard design criteria.

    Args:
        sequence: A DNA sequence (case-insensitive, whitespace ignored).

    Returns:
        A :class:`PrimerAnalysis` holding every measurement plus warnings for
        any criterion the primer fails.

    Raises:
        InvalidSequenceError: If ``sequence`` is not valid DNA.
    """
    seq = _require_valid_dna(sequence)

    length = len(seq)
    gc = gc_content(seq)
    tm = float(melting_temperature(seq))
    has_clamp = seq[-1] in _GC_BASES
    homopolymer = _longest_homopolymer_run(seq)
    self_comp = _self_complementarity(seq)

    warnings: list[str] = []

    if length < MIN_LENGTH:
        warnings.append(f"Too short ({length} bases; ideal {MIN_LENGTH}-{MAX_LENGTH}).")
    elif length > MAX_LENGTH:
        warnings.append(f"Too long ({length} bases; ideal {MIN_LENGTH}-{MAX_LENGTH}).")

    if gc < MIN_GC_PERCENT:
        warnings.append(
            f"GC content too low ({gc:.1f}%; ideal "
            f"{MIN_GC_PERCENT:.0f}-{MAX_GC_PERCENT:.0f}%)."
        )
    elif gc > MAX_GC_PERCENT:
        warnings.append(
            f"GC content too high ({gc:.1f}%; ideal "
            f"{MIN_GC_PERCENT:.0f}-{MAX_GC_PERCENT:.0f}%)."
        )

    if tm < MIN_MELTING_TEMP:
        warnings.append(
            f"Melting temperature too low ({tm:.1f} C; ideal "
            f"{MIN_MELTING_TEMP:.0f}-{MAX_MELTING_TEMP:.0f} C)."
        )
    elif tm > MAX_MELTING_TEMP:
        warnings.append(
            f"Melting temperature too high ({tm:.1f} C; ideal "
            f"{MIN_MELTING_TEMP:.0f}-{MAX_MELTING_TEMP:.0f} C)."
        )

    if not has_clamp:
        warnings.append("No GC clamp: the 3' end should be G or C.")

    if homopolymer > MAX_HOMOPOLYMER_RUN:
        warnings.append(
            f"Homopolymer run of {homopolymer} bases "
            f"(keep under {MAX_HOMOPOLYMER_RUN + 1})."
        )

    if self_comp > MAX_SELF_COMPLEMENTARITY:
        warnings.append(
            f"Self-complementary stretch of {self_comp} bases may form a "
            f"hairpin or primer-dimer."
        )

    return PrimerAnalysis(
        sequence=seq,
        length=length,
        gc_percent=gc,
        melting_temp=tm,
        has_gc_clamp=has_clamp,
        max_homopolymer_run=homopolymer,
        self_complementarity=self_comp,
        warnings=warnings,
    )
