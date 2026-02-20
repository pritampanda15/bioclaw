"""Start the FastAPI server."""

from __future__ import annotations

import click


@click.command()
@click.option("--host", "-h", default=None, help="Bind host")
@click.option("--port", "-p", default=None, type=int, help="Bind port")
@click.option("--reload", is_flag=True, help="Enable auto-reload for development")
def serve(host: str | None, port: int | None, reload: bool) -> None:
    """Start the BioClaw API server."""
    import uvicorn

    from bioclaw.config.settings import get_settings
    from bioclaw.utils.logging import setup_logging

    settings = get_settings()
    setup_logging(settings.log_level)

    uvicorn.run(
        "bioclaw.server.app:create_app",
        factory=True,
        host=host or settings.server.host,
        port=port or settings.server.port,
        reload=reload,
    )
