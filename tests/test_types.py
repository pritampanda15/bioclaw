"""Tests for tool type definitions."""

from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult


def test_tool_parameter_to_json_schema():
    param = ToolParameter(
        name="smiles",
        type=ParameterType.STRING,
        description="SMILES string",
        required=True,
    )
    schema = param.to_json_schema()
    assert schema == {"type": "string", "description": "SMILES string"}


def test_tool_parameter_with_enum():
    param = ToolParameter(
        name="mode",
        type=ParameterType.STRING,
        description="Mode",
        enum=["2d", "3d"],
    )
    schema = param.to_json_schema()
    assert schema["enum"] == ["2d", "3d"]


def test_tool_parameter_array_type():
    param = ToolParameter(
        name="smiles_list",
        type=ParameterType.ARRAY,
        description="List of SMILES",
        items_type=ParameterType.STRING,
    )
    schema = param.to_json_schema()
    assert schema["type"] == "array"
    assert schema["items"] == {"type": "string"}


def test_tool_metadata_to_claude_schema():
    meta = ToolMetadata(
        name="test_tool",
        description="A test tool",
        parameters=[
            ToolParameter(name="input", type=ParameterType.STRING, description="Input value"),
            ToolParameter(name="flag", type=ParameterType.BOOLEAN, description="A flag", required=False),
        ],
        category="test",
    )
    schema = meta.to_claude_tool_schema()
    assert schema["name"] == "test_tool"
    assert schema["description"] == "A test tool"
    assert "input" in schema["input_schema"]["properties"]
    assert "flag" in schema["input_schema"]["properties"]
    assert schema["input_schema"]["required"] == ["input"]


def test_tool_result_to_dict():
    result = ToolResult(
        success=True,
        data={"key": "value"},
        summary="Test result",
        warnings=["warning1"],
    )
    d = result.to_dict()
    assert d["success"] is True
    assert d["data"] == {"key": "value"}
    assert d["summary"] == "Test result"
    assert d["warnings"] == ["warning1"]


def test_tool_result_error():
    result = ToolResult(success=False, error="Something failed")
    d = result.to_dict()
    assert d["success"] is False
    assert d["error"] == "Something failed"
