"""Type definitions for the tool system."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any


class ParameterType(str, Enum):
    STRING = "string"
    INTEGER = "integer"
    NUMBER = "number"
    BOOLEAN = "boolean"
    ARRAY = "array"
    OBJECT = "object"


@dataclass
class ToolParameter:
    name: str
    type: ParameterType
    description: str
    required: bool = True
    default: Any = None
    enum: list[str] | None = None
    items_type: ParameterType | None = None  # For array types

    def to_json_schema(self) -> dict[str, Any]:
        schema: dict[str, Any] = {
            "type": self.type.value,
            "description": self.description,
        }
        if self.enum:
            schema["enum"] = self.enum
        if self.type == ParameterType.ARRAY and self.items_type:
            schema["items"] = {"type": self.items_type.value}
        if self.default is not None:
            schema["default"] = self.default
        return schema


@dataclass
class ToolMetadata:
    name: str
    description: str
    parameters: list[ToolParameter]
    category: str = "common"
    is_long_running: bool = False
    required_binaries: list[str] = field(default_factory=list)

    def to_claude_tool_schema(self) -> dict[str, Any]:
        """Convert to Claude API tool schema format."""
        properties = {}
        required = []
        for param in self.parameters:
            properties[param.name] = param.to_json_schema()
            if param.required:
                required.append(param.name)

        return {
            "name": self.name,
            "description": self.description,
            "input_schema": {
                "type": "object",
                "properties": properties,
                "required": required,
            },
        }


@dataclass
class ToolResult:
    success: bool
    data: dict[str, Any] = field(default_factory=dict)
    summary: str = ""
    files: list[Path] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "success": self.success,
            "data": self.data,
            "summary": self.summary,
        }
        if self.files:
            result["files"] = [str(f) for f in self.files]
        if self.warnings:
            result["warnings"] = self.warnings
        if self.error:
            result["error"] = self.error
        return result
