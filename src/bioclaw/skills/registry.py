"""Auto-discovery and registration of tools."""

from __future__ import annotations

import importlib
import pkgutil
from typing import Any

import bioclaw.skills
from bioclaw.skills.base import BaseTool
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.registry")


class ToolRegistry:
    """Discovers and manages all available tools."""

    def __init__(self) -> None:
        self._tools: dict[str, BaseTool] = {}

    def discover(self, categories: list[str] | None = None) -> None:
        """Auto-discover all BaseTool subclasses in skills subpackages."""
        package = bioclaw.skills
        for importer, modname, ispkg in pkgutil.walk_packages(
            package.__path__, prefix=package.__name__ + "."
        ):
            # Skip base/registry/types modules
            if modname.endswith((".base", ".registry", ".types")):
                continue

            # Filter by category if specified
            if categories:
                parts = modname.split(".")
                # e.g. bioclaw.skills.drugdiscovery.molecular_properties
                if len(parts) >= 4:
                    category = parts[-2]  # drugdiscovery, ngs, common
                    if category not in categories:
                        continue

            try:
                module = importlib.import_module(modname)
            except ImportError as e:
                logger.debug(f"Skipping {modname}: {e}")
                continue

            for attr_name in dir(module):
                attr = getattr(module, attr_name)
                if (
                    isinstance(attr, type)
                    and issubclass(attr, BaseTool)
                    and attr is not BaseTool
                ):
                    try:
                        instance = attr()
                        self._tools[instance.metadata.name] = instance
                        logger.debug(f"Registered tool: {instance.metadata.name}")
                    except Exception as e:
                        logger.debug(f"Failed to instantiate {attr_name}: {e}")

    def get(self, name: str) -> BaseTool | None:
        return self._tools.get(name)

    def list_tools(self) -> list[BaseTool]:
        return list(self._tools.values())

    def get_claude_tool_schemas(self) -> list[dict[str, Any]]:
        """Get all tool schemas in Claude API format."""
        return [tool.metadata.to_claude_tool_schema() for tool in self._tools.values()]

    def check_availability(self) -> dict[str, bool]:
        """Check which tools are available."""
        return {name: tool.check_available() for name, tool in self._tools.items()}


_registry: ToolRegistry | None = None


def get_registry(categories: list[str] | None = None) -> ToolRegistry:
    global _registry
    if _registry is None:
        _registry = ToolRegistry()
        _registry.discover(categories)
    return _registry
