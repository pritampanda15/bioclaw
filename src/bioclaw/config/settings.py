"""Application settings from env vars and YAML config."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml
from pydantic import Field
from pydantic_settings import BaseSettings


def _load_yaml_config() -> dict[str, Any]:
    """Load bioclaw.yaml from CWD or project root."""
    for candidate in [Path("bioclaw.yaml"), Path(__file__).parents[3] / "bioclaw.yaml"]:
        if candidate.exists():
            with open(candidate) as f:
                return yaml.safe_load(f) or {}
    return {}


class ServerSettings(BaseSettings):
    host: str = "127.0.0.1"
    port: int = 8000


class ToolSettings(BaseSettings):
    enabled_categories: list[str] = Field(default_factory=lambda: ["drugdiscovery", "ngs", "common"])


class Settings(BaseSettings):
    model_config = {"env_prefix": "BIOCLAW_"}

    anthropic_api_key: str = Field(default="", alias="ANTHROPIC_API_KEY")
    model: str = "claude-sonnet-4-5-20250929"
    max_tool_rounds: int = 15
    data_dir: Path = Path("./bioclaw_data")
    log_level: str = "INFO"
    server: ServerSettings = Field(default_factory=ServerSettings)
    tools: ToolSettings = Field(default_factory=ToolSettings)

    def model_post_init(self, __context: Any) -> None:
        # Overlay YAML config for nested fields
        yaml_cfg = _load_yaml_config()
        if "server" in yaml_cfg and isinstance(yaml_cfg["server"], dict):
            for k, v in yaml_cfg["server"].items():
                if hasattr(self.server, k):
                    object.__setattr__(self.server, k, v)
        if "tools" in yaml_cfg and isinstance(yaml_cfg["tools"], dict):
            for k, v in yaml_cfg["tools"].items():
                if hasattr(self.tools, k):
                    object.__setattr__(self.tools, k, v)

        self.data_dir.mkdir(parents=True, exist_ok=True)


_settings: Settings | None = None


def get_settings() -> Settings:
    global _settings
    if _settings is None:
        _settings = Settings()
    return _settings
