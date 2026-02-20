"""Job data models."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Any
from uuid import uuid4

from bioclaw.skills.types import ToolResult


class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class Job:
    tool_name: str
    parameters: dict[str, Any]
    id: str = field(default_factory=lambda: uuid4().hex)
    status: JobStatus = JobStatus.PENDING
    progress: int | None = None
    created_at: datetime = field(default_factory=datetime.now)
    started_at: datetime | None = None
    completed_at: datetime | None = None
    result: ToolResult | None = None
    error: str | None = None

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "id": self.id,
            "tool_name": self.tool_name,
            "parameters": self.parameters,
            "status": self.status.value,
            "progress": self.progress,
            "created_at": self.created_at.isoformat(),
        }
        if self.started_at:
            d["started_at"] = self.started_at.isoformat()
        if self.completed_at:
            d["completed_at"] = self.completed_at.isoformat()
        if self.result:
            d["result"] = self.result.to_dict()
        if self.error:
            d["error"] = self.error
        return d
