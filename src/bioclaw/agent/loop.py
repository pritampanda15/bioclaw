"""ReAct agent loop — reason, act, observe."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, AsyncGenerator

from bioclaw.agent.claude_client import ClaudeClient
from bioclaw.agent.memory import ConversationMemory
from bioclaw.agent.prompts import SYSTEM_PROMPT
from bioclaw.config.settings import get_settings
from bioclaw.skills.registry import ToolRegistry, get_registry
from bioclaw.utils.logging import get_logger

logger = get_logger("agent.loop")


@dataclass
class AgentEvent:
    """Event emitted during agent execution."""

    type: str  # text_delta, tool_call, tool_result, error, done
    data: Any = None


class AgentLoop:
    """ReAct agent loop using Claude's tool-use API."""

    def __init__(
        self,
        client: ClaudeClient | None = None,
        registry: ToolRegistry | None = None,
        memory: ConversationMemory | None = None,
    ) -> None:
        settings = get_settings()
        self.client = client or ClaudeClient()
        self.registry = registry or get_registry(settings.tools.enabled_categories)
        self.memory = memory or ConversationMemory()
        self.max_tool_rounds = settings.max_tool_rounds

    async def run(self, user_message: str) -> AsyncGenerator[AgentEvent, None]:
        """Process a user message through the agent loop.

        Yields AgentEvents as the agent reasons and acts.
        """
        self.memory.add_user_message(user_message)
        tools = self.registry.get_claude_tool_schemas()

        for round_num in range(self.max_tool_rounds):
            logger.debug(f"Agent round {round_num + 1}")

            try:
                response = self.client.create_message(
                    messages=self.memory.get_messages(),
                    system=SYSTEM_PROMPT,
                    tools=tools,
                )
            except Exception as e:
                yield AgentEvent(type="error", data=str(e))
                return

            # Collect content blocks
            content_blocks = response.content
            assistant_content = []
            tool_calls = []

            for block in content_blocks:
                if block.type == "text":
                    assistant_content.append({"type": "text", "text": block.text})
                    yield AgentEvent(type="text_delta", data=block.text)
                elif block.type == "tool_use":
                    assistant_content.append({
                        "type": "tool_use",
                        "id": block.id,
                        "name": block.name,
                        "input": block.input,
                    })
                    tool_calls.append(block)

            # Store assistant message
            self.memory.add_assistant_message(assistant_content)

            # If no tool calls, we're done
            if not tool_calls:
                yield AgentEvent(type="done")
                return

            # Execute tools and feed results back
            tool_results = []
            for tc in tool_calls:
                yield AgentEvent(
                    type="tool_call",
                    data={"name": tc.name, "input": tc.input, "id": tc.id},
                )

                result = await self._execute_tool(tc.name, tc.input)
                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": tc.id,
                    "content": json.dumps(result.to_dict()),
                })

                yield AgentEvent(
                    type="tool_result",
                    data={"name": tc.name, "result": result.to_dict()},
                )

            self.memory.add_tool_results(tool_results)

        # Hit max rounds
        yield AgentEvent(
            type="error",
            data=f"Reached maximum tool rounds ({self.max_tool_rounds})",
        )

    async def _execute_tool(self, name: str, inputs: dict[str, Any]):
        """Execute a tool by name."""
        from bioclaw.skills.types import ToolResult

        tool = self.registry.get(name)
        if tool is None:
            return ToolResult(
                success=False,
                error=f"Unknown tool: {name}",
                summary=f"Tool '{name}' not found in registry",
            )
        return await tool.safe_execute(**inputs)

    def reset(self) -> None:
        """Clear conversation memory."""
        self.memory.clear()
