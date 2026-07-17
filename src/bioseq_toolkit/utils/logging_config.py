"""Logging helpers for BioSeq Toolkit.

Design rule: libraries emit log records but do NOT configure global logging;
the application (our CLI) calls ``configure_logging()`` once at startup to
decide the level and format. Modules just call ``get_logger(__name__)``.
"""

from __future__ import annotations

import logging

_DEFAULT_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
_DEFAULT_DATEFMT = "%Y-%m-%d %H:%M:%S"


def get_logger(name: str) -> logging.Logger:
    """Return a module-level logger.

    Args:
        name: Usually ``__name__`` from the calling module, so log records
            are tagged with where they came from.

    Returns:
        A configured :class:`logging.Logger` instance.
    """
    return logging.getLogger(name)


def configure_logging(level: int = logging.INFO) -> None:
    """Configure root logging for the application.

    Call this ONCE, from the application entry point (the CLI), not from
    library modules. Safe to call repeatedly: it clears existing handlers
    first so logs are not duplicated.

    Args:
        level: The minimum severity to emit (e.g. ``logging.DEBUG``).
    """
    root = logging.getLogger()
    root.setLevel(level)

    # Remove old handlers so repeated calls don't stack up duplicate output.
    for handler in list(root.handlers):
        root.removeHandler(handler)

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT, _DEFAULT_DATEFMT))
    root.addHandler(handler)
