"""Main CLI entry point."""

from __future__ import annotations

import click

from bioclaw import __version__


@click.group()
@click.version_option(version=__version__, prog_name="bioclaw")
def cli() -> None:
    """BioClaw — Bioinformatics AI agent for drug discovery & NGS analysis."""
    pass


# Register sub-commands
from bioclaw.cli.chat import chat  # noqa: E402
from bioclaw.cli.jobs import jobs  # noqa: E402
from bioclaw.cli.run import run  # noqa: E402
from bioclaw.cli.serve import serve  # noqa: E402

cli.add_command(chat)
cli.add_command(run)
cli.add_command(jobs)
cli.add_command(serve)


@cli.group(name="tools")
def tools_group() -> None:
    """Manage available tools."""
    pass


@tools_group.command(name="list")
def tools_list() -> None:
    """List all registered tools."""
    from rich.console import Console
    from rich.table import Table

    from bioclaw.config.settings import get_settings
    from bioclaw.skills.registry import get_registry
    from bioclaw.utils.logging import setup_logging

    settings = get_settings()
    setup_logging(settings.log_level)
    registry = get_registry(settings.tools.enabled_categories)

    console = Console()
    table = Table(title="BioClaw Tools")
    table.add_column("Name", style="cyan")
    table.add_column("Category", style="green")
    table.add_column("Description")
    table.add_column("Available", style="bold")

    for tool in registry.list_tools():
        available = "✓" if tool.check_available() else "✗"
        style = "green" if tool.check_available() else "red"
        table.add_row(
            tool.metadata.name,
            tool.metadata.category,
            tool.metadata.description[:80] + "..." if len(tool.metadata.description) > 80 else tool.metadata.description,
            f"[{style}]{available}[/{style}]",
        )

    console.print(table)


@tools_group.command(name="check")
def tools_check() -> None:
    """Check tool availability and dependencies."""
    from rich.console import Console

    from bioclaw.config.settings import get_settings
    from bioclaw.skills.registry import get_registry
    from bioclaw.utils.logging import setup_logging

    settings = get_settings()
    setup_logging(settings.log_level)
    registry = get_registry(settings.tools.enabled_categories)

    console = Console()
    availability = registry.check_availability()

    available = sum(1 for v in availability.values() if v)
    total = len(availability)
    console.print(f"\n[bold]Tool Availability: {available}/{total}[/bold]\n")

    for name, is_available in sorted(availability.items()):
        tool = registry.get(name)
        if is_available:
            console.print(f"  [green]✓[/green] {name}")
        else:
            missing = tool.metadata.required_binaries if tool else []
            console.print(f"  [red]✗[/red] {name} (missing: {', '.join(missing) if missing else 'library'})")
