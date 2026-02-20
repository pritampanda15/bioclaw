"""Job management commands."""

from __future__ import annotations

import asyncio

import click
from rich.console import Console
from rich.table import Table


@click.group()
def jobs() -> None:
    """Manage background jobs."""
    pass


@jobs.command(name="list")
@click.option("--status", "-s", type=click.Choice(["all", "running", "completed", "failed"]), default="all")
def jobs_list(status: str) -> None:
    """List all jobs."""
    console = Console()

    async def _list():
        from bioclaw.jobs.manager import get_job_manager

        manager = get_job_manager()
        all_jobs = await manager.list_jobs(status_filter=None if status == "all" else status)

        if not all_jobs:
            console.print("[dim]No jobs found.[/dim]")
            return

        table = Table(title="Jobs")
        table.add_column("ID", style="cyan")
        table.add_column("Tool", style="green")
        table.add_column("Status")
        table.add_column("Created")
        table.add_column("Progress")

        for job in all_jobs:
            status_style = {
                "running": "yellow",
                "completed": "green",
                "failed": "red",
                "cancelled": "dim",
                "pending": "blue",
            }.get(job.status, "white")
            table.add_row(
                job.id[:8],
                job.tool_name,
                f"[{status_style}]{job.status}[/{status_style}]",
                job.created_at.strftime("%H:%M:%S"),
                f"{job.progress}%" if job.progress is not None else "-",
            )

        console.print(table)

    asyncio.run(_list())


@jobs.command(name="status")
@click.argument("job_id")
def jobs_status(job_id: str) -> None:
    """Check status of a specific job."""
    console = Console()

    async def _status():
        from bioclaw.jobs.manager import get_job_manager

        manager = get_job_manager()
        job = await manager.get_job(job_id)
        if job is None:
            console.print(f"[red]Job not found:[/red] {job_id}")
            return
        console.print(f"[bold]Job {job.id[:8]}[/bold]")
        console.print(f"  Tool: {job.tool_name}")
        console.print(f"  Status: {job.status}")
        console.print(f"  Created: {job.created_at}")
        if job.result:
            console.print(f"  Result: {job.result.summary}")

    asyncio.run(_status())


@jobs.command(name="cancel")
@click.argument("job_id")
def jobs_cancel(job_id: str) -> None:
    """Cancel a running job."""
    console = Console()

    async def _cancel():
        from bioclaw.jobs.manager import get_job_manager

        manager = get_job_manager()
        success = await manager.cancel_job(job_id)
        if success:
            console.print(f"[green]Job {job_id[:8]} cancelled.[/green]")
        else:
            console.print(f"[red]Failed to cancel job {job_id[:8]}.[/red]")

    asyncio.run(_cancel())
