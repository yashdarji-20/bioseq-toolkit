"""Custom exceptions for BioSeq Toolkit.

Defining domain-specific exceptions (rather than raising bare ``ValueError``)
lets callers catch exactly the failure they care about:

    try:
        gc = gc_content(user_input)
    except InvalidSequenceError:
        ...  # handle only *our* sequence problem, nothing else

All toolkit exceptions inherit from :class:`BioSeqError`, so a caller can also
catch every error this library raises with a single ``except BioSeqError``.
"""

from __future__ import annotations


class BioSeqError(Exception):
    """Base class for all BioSeq Toolkit errors."""


class InvalidSequenceError(BioSeqError):
    """Raised when a sequence contains characters outside its alphabet."""