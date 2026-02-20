"""Base class for all BioClaw tools."""

from __future__ import annotations

import shutil
from abc import ABC, abstractmethod

from bioclaw.skills.types import ToolMetadata, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.base")


class BaseTool(ABC):
    """Abstract base class for all tools."""

    @property
    @abstractmethod
    def metadata(self) -> ToolMetadata:
        """Return tool metadata including name, description, and parameters."""
        ...

    @abstractmethod
    async def execute(self, **kwargs) -> ToolResult:
        """Execute the tool with given parameters."""
        ...

    def check_available(self) -> bool:
        """Check if required binaries/libraries are available."""
        for binary in self.metadata.required_binaries:
            if shutil.which(binary) is None:
                logger.debug(f"Binary not found: {binary}")
                return False
        return True

    async def safe_execute(self, **kwargs) -> ToolResult:
        """Execute with error handling."""
        try:
            if not self.check_available():
                missing = [b for b in self.metadata.required_binaries if shutil.which(b) is None]
                return ToolResult(
                    success=False,
                    error=f"Missing required binaries: {', '.join(missing)}",
                    summary=f"Tool unavailable: missing {', '.join(missing)}",
                )
            return await self.execute(**kwargs)
        except Exception as e:
            logger.exception(f"Tool {self.metadata.name} failed")
            return ToolResult(
                success=False,
                error=str(e),
                summary=f"Tool {self.metadata.name} failed: {e}",
            )
