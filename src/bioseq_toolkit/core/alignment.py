"""Global sequence alignment (Needleman-Wunsch) for BioSeq Toolkit.

Implements the classic dynamic-programming algorithm for *global* alignment:
the two sequences are aligned end to end, inserting gaps where needed to
maximise a simple match / mismatch / gap score.

The algorithm builds an ``(m+1) x (n+1)`` scoring matrix where cell ``[i][j]``
holds the best achievable score for aligning the first ``i`` bases of one
sequence with the first ``j`` of the other. Each cell depends on only three
neighbours, which is what turns an exponential search into ``O(m*n)`` work::

    score[i][j] = max(
        score[i-1][j-1] + match_or_mismatch,   # align both bases
        score[i-1][j]   + gap,                 # gap in sequence 2
        score[i][j-1]   + gap,                 # gap in sequence 1
    )

A traceback from the bottom-right corner then reconstructs the alignment.
"""

from __future__ import annotations

from dataclasses import dataclass

# Default scoring scheme. Gaps cost more than mismatches because insertions
# and deletions are rarer evolutionary events than substitutions.
_DEFAULT_MATCH = 1
_DEFAULT_MISMATCH = -1
_DEFAULT_GAP = -2

_GAP_CHARACTER = "-"


@dataclass(frozen=True)
class Alignment:
    """The result of aligning two sequences.

    Attributes:
        aligned_seq1: First sequence with gap characters inserted.
        aligned_seq2: Second sequence with gap characters inserted.
        score: Total alignment score under the chosen scoring scheme.
    """

    aligned_seq1: str
    aligned_seq2: str
    score: int

    @property
    def identity(self) -> float:
        """Percentage of aligned columns where both bases are identical."""
        if not self.aligned_seq1:
            return 0.0
        matches = sum(
            1
            for a, b in zip(self.aligned_seq1, self.aligned_seq2)
            if a == b and a != _GAP_CHARACTER
        )
        return matches / len(self.aligned_seq1) * 100

    def formatted(self) -> str:
        """Return a three-line human-readable alignment view."""
        middle = "".join(
            "|" if a == b and a != _GAP_CHARACTER else " "
            for a, b in zip(self.aligned_seq1, self.aligned_seq2)
        )
        return f"{self.aligned_seq1}\n{middle}\n{self.aligned_seq2}"


def _build_matrix(
    seq1: str, seq2: str, match: int, mismatch: int, gap: int
) -> list[list[int]]:
    """Build the Needleman-Wunsch scoring matrix.

    Args:
        seq1: First sequence (indexes the rows).
        seq2: Second sequence (indexes the columns).
        match: Score awarded when two bases are identical.
        mismatch: Score applied when two bases differ.
        gap: Penalty applied for a gap.

    Returns:
        An ``(len(seq1)+1) x (len(seq2)+1)`` matrix of scores.
    """
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = [[0] * cols for _ in range(rows)]

    # First column/row: aligning a prefix against nothing is all gaps.
    for i in range(1, rows):
        matrix[i][0] = matrix[i - 1][0] + gap
    for j in range(1, cols):
        matrix[0][j] = matrix[0][j - 1] + gap

    for i in range(1, rows):
        for j in range(1, cols):
            pair_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            matrix[i][j] = max(
                matrix[i - 1][j - 1] + pair_score,  # align both bases
                matrix[i - 1][j] + gap,  # gap in seq2
                matrix[i][j - 1] + gap,  # gap in seq1
            )

    return matrix


def _traceback(
    matrix: list[list[int]],
    seq1: str,
    seq2: str,
    match: int,
    mismatch: int,
    gap: int,
) -> tuple[str, str]:
    """Reconstruct the optimal alignment by walking back through the matrix.

    Starts at the bottom-right corner and at each step determines which
    neighbour produced the current score, moving there and recording the
    corresponding alignment column.

    Returns:
        A tuple of ``(aligned_seq1, aligned_seq2)``.
    """
    aligned1: list[str] = []
    aligned2: list[str] = []
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        pair_score = (
            match if i > 0 and j > 0 and seq1[i - 1] == seq2[j - 1] else mismatch
        )
        if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + pair_score:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap:
            aligned1.append(seq1[i - 1])
            aligned2.append(_GAP_CHARACTER)
            i -= 1
        else:
            aligned1.append(_GAP_CHARACTER)
            aligned2.append(seq2[j - 1])
            j -= 1

    # The traceback builds the alignment backwards, so reverse it.
    return "".join(reversed(aligned1)), "".join(reversed(aligned2))


def align(
    seq1: str,
    seq2: str,
    *,
    match: int = _DEFAULT_MATCH,
    mismatch: int = _DEFAULT_MISMATCH,
    gap: int = _DEFAULT_GAP,
) -> Alignment:
    """Globally align two sequences using the Needleman-Wunsch algorithm.

    Args:
        seq1: First sequence.
        seq2: Second sequence.
        match: Score for identical bases (default ``1``).
        mismatch: Score for differing bases (default ``-1``).
        gap: Penalty for a gap (default ``-2``).

    Returns:
        An :class:`Alignment` holding both gapped sequences and the score.
    """
    s1 = seq1.strip().upper()
    s2 = seq2.strip().upper()

    matrix = _build_matrix(s1, s2, match, mismatch, gap)
    aligned1, aligned2 = _traceback(matrix, s1, s2, match, mismatch, gap)

    return Alignment(
        aligned_seq1=aligned1,
        aligned_seq2=aligned2,
        score=matrix[len(s1)][len(s2)],
    )
