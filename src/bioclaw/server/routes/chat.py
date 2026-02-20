"""Chat API routes and WebSocket handler."""

from __future__ import annotations

from fastapi import APIRouter, WebSocket, WebSocketDisconnect
from pydantic import BaseModel

from bioclaw.agent.loop import AgentLoop
from bioclaw.server.ws import ws_manager
from bioclaw.utils.logging import get_logger

logger = get_logger("server.routes.chat")
router = APIRouter(tags=["chat"])

# Per-session agent loops (simple in-memory for now)
_sessions: dict[str, AgentLoop] = {}


def _get_agent(session_id: str = "default") -> AgentLoop:
    if session_id not in _sessions:
        _sessions[session_id] = AgentLoop()
    return _sessions[session_id]


class ChatRequest(BaseModel):
    message: str
    session_id: str = "default"


class ChatResponse(BaseModel):
    response: str
    tool_calls: list[dict] = []
    tool_results: list[dict] = []


@router.post("/chat", response_model=ChatResponse)
async def chat(request: ChatRequest) -> ChatResponse:
    """Send a message and get a complete response."""
    agent = _get_agent(request.session_id)
    response_text = ""
    tool_calls = []
    tool_results = []

    async for event in agent.run(request.message):
        if event.type == "text_delta":
            response_text += event.data
        elif event.type == "tool_call":
            tool_calls.append(event.data)
        elif event.type == "tool_result":
            tool_results.append(event.data)
        elif event.type == "error":
            response_text += f"\n\nError: {event.data}"

    return ChatResponse(
        response=response_text,
        tool_calls=tool_calls,
        tool_results=tool_results,
    )


async def websocket_chat(websocket: WebSocket) -> None:
    """WebSocket handler for streaming chat."""
    await ws_manager.connect(websocket)
    session_id = "ws_default"

    try:
        while True:
            data = await websocket.receive_json()
            message = data.get("message", "")
            session_id = data.get("session_id", session_id)

            if data.get("action") == "clear":
                agent = _get_agent(session_id)
                agent.reset()
                await ws_manager.send_event(websocket, "cleared")
                continue

            agent = _get_agent(session_id)

            async for event in agent.run(message):
                await ws_manager.send_event(websocket, event.type, event.data)

    except WebSocketDisconnect:
        ws_manager.disconnect(websocket)
    except Exception as e:
        logger.exception("WebSocket error")
        try:
            await ws_manager.send_event(websocket, "error", str(e))
        except Exception:
            pass
        ws_manager.disconnect(websocket)
