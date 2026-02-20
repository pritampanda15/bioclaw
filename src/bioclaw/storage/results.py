"""Result file management."""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path
from uuid import uuid4

from bioclaw.config.settings import get_settings
from bioclaw.utils.logging import get_logger

logger = get_logger("storage.results")


class ResultManager:
    """Manages result files and metadata."""

    def __init__(self) -> None:
        settings = get_settings()
        self.results_dir = settings.data_dir / "results"
        self.results_dir.mkdir(parents=True, exist_ok=True)

    def create_result_dir(self, tool_name: str) -> Path:
        """Create a unique directory for tool results."""
        result_id = uuid4().hex[:12]
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        dir_name = f"{timestamp}_{tool_name}_{result_id}"
        result_dir = self.results_dir / dir_name
        result_dir.mkdir(parents=True, exist_ok=True)
        return result_dir

    def save_result(self, result_dir: Path, data: dict, filename: str = "result.json") -> Path:
        """Save result data as JSON."""
        output = result_dir / filename
        with open(output, "w") as f:
            json.dump(data, f, indent=2, default=str)
        return output

    def list_results(self, tool_name: str | None = None) -> list[Path]:
        """List all result directories, optionally filtered by tool name."""
        dirs = sorted(self.results_dir.iterdir(), reverse=True)
        if tool_name:
            dirs = [d for d in dirs if tool_name in d.name]
        return dirs

    def get_result(self, result_id: str) -> dict | None:
        """Load a result by ID (prefix match)."""
        for d in self.results_dir.iterdir():
            if result_id in d.name:
                result_file = d / "result.json"
                if result_file.exists():
                    with open(result_file) as f:
                        return json.load(f)
        return None

    def cleanup_old(self, max_age_days: int = 30) -> int:
        """Remove results older than max_age_days."""
        cutoff = datetime.now().timestamp() - (max_age_days * 86400)
        removed = 0
        for d in self.results_dir.iterdir():
            if d.is_dir() and d.stat().st_mtime < cutoff:
                shutil.rmtree(d)
                removed += 1
        return removed
