"""FastAPI application factory."""

from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from bioclaw import __version__
from bioclaw.config.settings import get_settings
from bioclaw.utils.logging import setup_logging


def create_app() -> FastAPI:
    settings = get_settings()
    setup_logging(settings.log_level)

    app = FastAPI(
        title="BioClaw API",
        version=__version__,
        description="Bioinformatics AI agent API",
    )

    app.add_middleware(
        CORSMiddleware,
        allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Register routes
    from bioclaw.server.routes.chat import router as chat_router
    from bioclaw.server.routes.files import router as files_router
    from bioclaw.server.routes.jobs import router as jobs_router
    from bioclaw.server.routes.tools import router as tools_router

    app.include_router(chat_router, prefix="/api")
    app.include_router(tools_router, prefix="/api")
    app.include_router(jobs_router, prefix="/api")
    app.include_router(files_router, prefix="/api")

    # WebSocket route
    from bioclaw.server.routes.chat import websocket_chat

    app.add_api_websocket_route("/ws/chat", websocket_chat)

    @app.get("/api/health")
    async def health():
        return {"status": "ok", "version": __version__}

    return app
