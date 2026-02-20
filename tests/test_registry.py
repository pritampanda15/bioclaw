"""Tests for tool registry."""

from bioclaw.skills.base import BaseTool
from bioclaw.skills.registry import ToolRegistry
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult


class DummyTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="dummy_tool",
            description="A dummy tool for testing",
            parameters=[
                ToolParameter(name="input", type=ParameterType.STRING, description="Test input"),
            ],
            category="test",
        )

    async def execute(self, **kwargs) -> ToolResult:
        return ToolResult(
            success=True,
            data={"echo": kwargs.get("input", "")},
            summary=f"Echo: {kwargs.get('input', '')}",
        )


def test_registry_manual_registration():
    registry = ToolRegistry()
    tool = DummyTool()
    registry._tools[tool.metadata.name] = tool

    assert registry.get("dummy_tool") is not None
    assert registry.get("nonexistent") is None
    assert len(registry.list_tools()) == 1


def test_registry_claude_schemas():
    registry = ToolRegistry()
    tool = DummyTool()
    registry._tools[tool.metadata.name] = tool

    schemas = registry.get_claude_tool_schemas()
    assert len(schemas) == 1
    assert schemas[0]["name"] == "dummy_tool"


def test_registry_check_availability():
    registry = ToolRegistry()
    tool = DummyTool()
    registry._tools[tool.metadata.name] = tool

    avail = registry.check_availability()
    assert avail["dummy_tool"] is True


def test_registry_discover():
    """Test auto-discovery finds at least one tool."""
    registry = ToolRegistry()
    registry.discover()
    # Should find at least file_converter and database_query (no external deps)
    tools = registry.list_tools()
    names = [t.metadata.name for t in tools]
    assert "file_converter" in names
    assert "database_query" in names
