"""Database query tool for searching local results."""

from __future__ import annotations

import json
from pathlib import Path

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult


class DatabaseQueryTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="database_query",
            description=(
                "Search and query local BioClaw results database. "
                "Can list recent results, search by tool name, or retrieve specific results."
            ),
            parameters=[
                ToolParameter(
                    name="action",
                    type=ParameterType.STRING,
                    description="Query action to perform",
                    enum=["list_recent", "search", "get"],
                ),
                ToolParameter(
                    name="tool_name",
                    type=ParameterType.STRING,
                    description="Filter by tool name (for list_recent/search)",
                    required=False,
                ),
                ToolParameter(
                    name="result_id",
                    type=ParameterType.STRING,
                    description="Result ID to retrieve (for get action)",
                    required=False,
                ),
                ToolParameter(
                    name="limit",
                    type=ParameterType.INTEGER,
                    description="Maximum number of results to return",
                    required=False,
                    default=10,
                ),
            ],
            category="common",
        )

    async def execute(self, **kwargs) -> ToolResult:
        from bioclaw.storage.results import ResultManager

        action = kwargs["action"]
        manager = ResultManager()

        if action == "list_recent":
            tool_name = kwargs.get("tool_name")
            limit = kwargs.get("limit", 10)
            dirs = manager.list_results(tool_name)[:limit]
            results = [{"name": d.name, "path": str(d)} for d in dirs]
            return ToolResult(
                success=True,
                data={"results": results, "count": len(results)},
                summary=f"Found {len(results)} recent results"
                + (f" for {tool_name}" if tool_name else ""),
            )

        elif action == "get":
            result_id = kwargs.get("result_id")
            if not result_id:
                return ToolResult(success=False, error="result_id required for get action")
            data = manager.get_result(result_id)
            if data is None:
                return ToolResult(success=False, error=f"Result not found: {result_id}")
            return ToolResult(
                success=True,
                data=data,
                summary=f"Retrieved result {result_id}",
            )

        elif action == "search":
            tool_name = kwargs.get("tool_name")
            dirs = manager.list_results(tool_name)
            results = []
            for d in dirs[:kwargs.get("limit", 10)]:
                result_file = d / "result.json"
                if result_file.exists():
                    with open(result_file) as f:
                        try:
                            result_data = json.load(f)
                            results.append({"name": d.name, "summary": result_data.get("summary", "")})
                        except json.JSONDecodeError:
                            pass
            return ToolResult(
                success=True,
                data={"results": results, "count": len(results)},
                summary=f"Found {len(results)} results matching query",
            )

        return ToolResult(success=False, error=f"Unknown action: {action}")
