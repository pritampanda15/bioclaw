"""Conversation memory management."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class ConversationMemory:
    """Manages conversation history for the agent loop."""

    messages: list[dict[str, Any]] = field(default_factory=list)
    max_messages: int = 100

    def add_user_message(self, content: str) -> None:
        self.messages.append({"role": "user", "content": content})
        self._trim()

    def add_assistant_message(self, content: list[dict[str, Any]]) -> None:
        """Add assistant message with content blocks (text + tool_use)."""
        self.messages.append({"role": "assistant", "content": content})
        self._trim()

    def add_tool_results(self, results: list[dict[str, Any]]) -> None:
        """Add tool result blocks as a user message."""
        self.messages.append({"role": "user", "content": results})
        self._trim()

    def get_messages(self) -> list[dict[str, Any]]:
        return list(self.messages)

    def clear(self) -> None:
        self.messages.clear()

    def _trim(self) -> None:
        if len(self.messages) > self.max_messages:
            # Keep first message (context) and trim oldest after that
            self.messages = self.messages[-self.max_messages :]
