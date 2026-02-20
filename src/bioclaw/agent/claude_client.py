"""Anthropic API wrapper for BioClaw."""

from __future__ import annotations

from typing import Any

import anthropic

from bioclaw.config.settings import get_settings
from bioclaw.utils.logging import get_logger

logger = get_logger("agent.claude_client")


class ClaudeClient:
    """Wrapper around the Anthropic API."""

    def __init__(self, api_key: str | None = None, model: str | None = None) -> None:
        settings = get_settings()
        self.api_key = api_key or settings.anthropic_api_key
        self.model = model or settings.model
        self.client = anthropic.Anthropic(api_key=self.api_key)

    def create_message(
        self,
        messages: list[dict[str, Any]],
        system: str,
        tools: list[dict[str, Any]],
        max_tokens: int = 4096,
    ) -> anthropic.types.Message:
        """Create a non-streaming message."""
        return self.client.messages.create(
            model=self.model,
            max_tokens=max_tokens,
            system=system,
            messages=messages,
            tools=tools if tools else anthropic.NOT_GIVEN,
        )

    def create_message_stream(
        self,
        messages: list[dict[str, Any]],
        system: str,
        tools: list[dict[str, Any]],
        max_tokens: int = 4096,
    ):
        """Create a streaming message, returning a context manager."""
        return self.client.messages.stream(
            model=self.model,
            max_tokens=max_tokens,
            system=system,
            messages=messages,
            tools=tools if tools else anthropic.NOT_GIVEN,
        )
