"""Tool listing API routes."""

from __future__ import annotations

from fastapi import APIRouter

from bioclaw.config.settings import get_settings
from bioclaw.skills.registry import get_registry

router = APIRouter(tags=["tools"])


@router.get("/tools")
async def list_tools() -> list[dict]:
    """List all available tools with their metadata."""
    settings = get_settings()
    registry = get_registry(settings.tools.enabled_categories)
    tools = []
    for tool in registry.list_tools():
        meta = tool.metadata
        tools.append({
            "name": meta.name,
            "description": meta.description,
            "category": meta.category,
            "is_long_running": meta.is_long_running,
            "available": tool.check_available(),
            "parameters": [
                {
                    "name": p.name,
                    "type": p.type.value,
                    "description": p.description,
                    "required": p.required,
                }
                for p in meta.parameters
            ],
        })
    return tools


@router.get("/tools/check")
async def check_tools() -> dict:
    """Check availability of all tools."""
    settings = get_settings()
    registry = get_registry(settings.tools.enabled_categories)
    availability = registry.check_availability()
    return {
        "total": len(availability),
        "available": sum(1 for v in availability.values() if v),
        "tools": availability,
    }
