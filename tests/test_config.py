"""Tests for configuration."""

from pathlib import Path

from bioclaw.config.settings import Settings, get_settings


def test_settings_defaults():
    settings = get_settings()
    assert settings.model == "claude-sonnet-4-5-20250929"
    assert settings.max_tool_rounds == 15
    assert settings.log_level == "DEBUG"  # Set in conftest


def test_settings_data_dir_created(tmp_path: Path):
    settings = get_settings()
    assert settings.data_dir.exists()


def test_settings_server_defaults():
    settings = get_settings()
    assert settings.server.host == "127.0.0.1"
    assert settings.server.port == 8000
