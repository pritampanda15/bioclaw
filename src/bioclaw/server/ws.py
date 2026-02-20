"""WebSocket connection manager."""

from __future__ import annotations

import json
from typing import Any

from fastapi import WebSocket

from bioclaw.utils.logging import get_logger

logger = get_logger("server.ws")


class ConnectionManager:
    """Manages WebSocket connections."""

    def __init__(self) -> None:
        self.active_connections: list[WebSocket] = []

    async def connect(self, websocket: WebSocket) -> None:
        await websocket.accept()
        self.active_connections.append(websocket)
        logger.info(f"WebSocket connected. Total: {len(self.active_connections)}")

    def disconnect(self, websocket: WebSocket) -> None:
        self.active_connections.remove(websocket)
        logger.info(f"WebSocket disconnected. Total: {len(self.active_connections)}")

    async def send_event(self, websocket: WebSocket, event_type: str, data: Any = None) -> None:
        """Send a typed event to a single client."""
        await websocket.send_json({"type": event_type, "data": data})

    async def broadcast(self, event_type: str, data: Any = None) -> None:
        """Broadcast an event to all connected clients."""
        message = json.dumps({"type": event_type, "data": data})
        for connection in self.active_connections:
            try:
                await connection.send_text(message)
            except Exception:
                pass


ws_manager = ConnectionManager()
