"""Input validators for BioClaw."""

from __future__ import annotations

import re


def validate_smiles(smiles: str) -> bool:
    """Basic SMILES syntax validation. For full validation, use RDKit."""
    if not smiles or not isinstance(smiles, str):
        return False
    # Basic check: only valid SMILES characters
    return bool(re.match(r"^[A-Za-z0-9@+\-\[\]()\\/%=#$.:~*]+$", smiles))


def validate_file_path(path: str, allowed_extensions: list[str] | None = None) -> bool:
    """Validate a file path string."""
    if not path or ".." in path:
        return False
    if allowed_extensions:
        return any(path.lower().endswith(ext) for ext in allowed_extensions)
    return True
