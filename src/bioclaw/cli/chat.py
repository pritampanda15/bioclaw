"""Interactive chat REPL."""

from __future__ import annotations

import asyncio

import click
from prompt_toolkit import PromptSession
from prompt_toolkit.history import InMemoryHistory
from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel

from bioclaw.agent.loop import AgentLoop
from bioclaw.config.settings import get_settings
from bioclaw.utils.logging import setup_logging


async def _chat_loop() -> None:
    settings = get_settings()
    setup_logging(settings.log_level)
    console = Console()

    console.print(
        Panel(
            "[bold cyan]BioClaw[/bold cyan] — Bioinformatics AI Agent\n"
            "Type your question or 'quit' to exit.",
            title="Welcome",
        )
    )

    agent = AgentLoop()
    session: PromptSession = PromptSession(history=InMemoryHistory())

    while True:
        try:
            user_input = await asyncio.get_event_loop().run_in_executor(
                None, lambda: session.prompt("\n🧬 You: ")
            )
        except (EOFError, KeyboardInterrupt):
            console.print("\n[dim]Goodbye![/dim]")
            break

        user_input = user_input.strip()
        if not user_input:
            continue
        if user_input.lower() in ("quit", "exit", "q"):
            console.print("[dim]Goodbye![/dim]")
            break
        if user_input.lower() == "/clear":
            agent.reset()
            console.print("[dim]Conversation cleared.[/dim]")
            continue

        console.print()
        text_buffer = ""

        async for event in agent.run(user_input):
            if event.type == "text_delta":
                text_buffer += event.data
            elif event.type == "tool_call":
                if text_buffer:
                    console.print(Markdown(text_buffer))
                    text_buffer = ""
                console.print(
                    f"  [yellow]⚙ Calling:[/yellow] [bold]{event.data['name']}[/bold]"
                    f"({_format_args(event.data['input'])})"
                )
            elif event.type == "tool_result":
                result = event.data["result"]
                if result["success"]:
                    console.print(f"  [green]✓[/green] {result['summary'][:200]}")
                else:
                    console.print(f"  [red]✗[/red] {result.get('error', 'Unknown error')}")
            elif event.type == "error":
                console.print(f"[red]Error:[/red] {event.data}")
            elif event.type == "done":
                if text_buffer:
                    console.print(Markdown(text_buffer))
                    text_buffer = ""

        # Flush any remaining text
        if text_buffer:
            console.print(Markdown(text_buffer))


def _format_args(args: dict) -> str:
    """Format tool arguments for display."""
    parts = []
    for k, v in args.items():
        val = repr(v) if isinstance(v, str) else str(v)
        if len(val) > 50:
            val = val[:47] + "..."
        parts.append(f"{k}={val}")
    return ", ".join(parts)


@click.command()
def chat() -> None:
    """Start interactive chat with BioClaw agent."""
    asyncio.run(_chat_loop())
