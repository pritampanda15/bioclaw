"""Logging setup for BioClaw."""

from __future__ import annotations

import logging
import sys

from rich.logging import RichHandler


def setup_logging(level: str = "INFO") -> None:
    """Configure root logger with Rich handler."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
    )
    # Quiet noisy libraries
    for name in ("httpx", "httpcore", "anthropic"):
        logging.getLogger(name).setLevel(logging.WARNING)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(f"bioclaw.{name}")
