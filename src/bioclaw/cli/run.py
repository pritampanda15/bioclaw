"""Direct tool execution command."""

from __future__ import annotations

import asyncio
import json

import click
from rich.console import Console
from rich.syntax import Syntax

from bioclaw.config.settings import get_settings
from bioclaw.skills.registry import get_registry
from bioclaw.utils.logging import setup_logging


@click.command()
@click.argument("tool_name")
@click.option("--params", "-p", type=str, help="JSON string of tool parameters")
@click.option("--param", "-P", type=(str, str), multiple=True, help="Key-value parameter pairs")
def run(tool_name: str, params: str | None, param: tuple[tuple[str, str], ...]) -> None:
    """Run a tool directly by name.

    Example: bioclaw run molecular_properties -P smiles "CC(=O)Oc1ccccc1C(=O)O"
    """
    settings = get_settings()
    setup_logging(settings.log_level)
    console = Console()

    registry = get_registry(settings.tools.enabled_categories)
    tool = registry.get(tool_name)

    if tool is None:
        console.print(f"[red]Tool not found:[/red] {tool_name}")
        console.print("Run 'bioclaw tools list' to see available tools.")
        raise SystemExit(1)

    # Build kwargs
    kwargs: dict = {}
    if params:
        kwargs = json.loads(params)
    for key, value in param:
        # Try to parse as JSON value, fall back to string
        try:
            kwargs[key] = json.loads(value)
        except json.JSONDecodeError:
            kwargs[key] = value

    async def _run():
        result = await tool.safe_execute(**kwargs)
        if result.success:
            console.print(f"[green]✓[/green] {result.summary}\n")
            console.print(Syntax(json.dumps(result.data, indent=2), "json"))
            if result.files:
                console.print(f"\n[bold]Output files:[/bold]")
                for f in result.files:
                    console.print(f"  {f}")
        else:
            console.print(f"[red]✗ Error:[/red] {result.error}")

        if result.warnings:
            for w in result.warnings:
                console.print(f"[yellow]⚠ {w}[/yellow]")

    asyncio.run(_run())
