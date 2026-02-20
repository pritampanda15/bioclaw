"""Shared test fixtures."""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock

import pytest


@pytest.fixture(autouse=True)
def _set_test_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Set up test environment variables."""
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant-test-key")
    monkeypatch.setenv("BIOCLAW_DATA_DIR", str(tmp_path / "bioclaw_data"))
    monkeypatch.setenv("BIOCLAW_LOG_LEVEL", "DEBUG")
    # Reset singleton settings
    import bioclaw.config.settings as settings_mod
    settings_mod._settings = None
    import bioclaw.skills.registry as registry_mod
    registry_mod._registry = None


@pytest.fixture
def sample_smiles() -> dict[str, str]:
    """Common molecule SMILES for testing."""
    return {
        "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
        "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "ethanol": "CCO",
        "benzene": "c1ccccc1",
    }
