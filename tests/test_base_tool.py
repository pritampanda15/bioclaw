"""Tests for base tool functionality."""

import pytest

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult


class FailingTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="failing_tool",
            description="A tool that always fails",
            parameters=[],
            category="test",
        )

    async def execute(self, **kwargs) -> ToolResult:
        raise ValueError("Intentional failure")


class MissingBinaryTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="missing_binary_tool",
            description="A tool with missing binary",
            parameters=[],
            category="test",
            required_binaries=["nonexistent_binary_xyz"],
        )

    async def execute(self, **kwargs) -> ToolResult:
        return ToolResult(success=True, summary="Should not reach here")


@pytest.mark.asyncio
async def test_safe_execute_catches_exceptions():
    tool = FailingTool()
    result = await tool.safe_execute()
    assert result.success is False
    assert "Intentional failure" in result.error


@pytest.mark.asyncio
async def test_safe_execute_checks_binaries():
    tool = MissingBinaryTool()
    result = await tool.safe_execute()
    assert result.success is False
    assert "nonexistent_binary_xyz" in result.error


def test_check_available_no_binaries():
    tool = FailingTool()
    assert tool.check_available() is True


def test_check_available_missing_binary():
    tool = MissingBinaryTool()
    assert tool.check_available() is False
