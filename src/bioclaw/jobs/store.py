"""SQLite persistence for jobs."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

import aiosqlite

from bioclaw.jobs.models import Job, JobStatus
from bioclaw.skills.types import ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("jobs.store")

CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS jobs (
    id TEXT PRIMARY KEY,
    tool_name TEXT NOT NULL,
    parameters TEXT NOT NULL,
    status TEXT NOT NULL DEFAULT 'pending',
    progress INTEGER,
    created_at TEXT NOT NULL,
    started_at TEXT,
    completed_at TEXT,
    result TEXT,
    error TEXT
)
"""


class JobStore:
    """SQLite-backed job persistence."""

    def __init__(self, db_path: Path) -> None:
        self.db_path = db_path
        self._initialized = False

    async def _ensure_db(self) -> None:
        if not self._initialized:
            self.db_path.parent.mkdir(parents=True, exist_ok=True)
            async with aiosqlite.connect(self.db_path) as db:
                await db.execute(CREATE_TABLE)
                await db.commit()
            self._initialized = True

    async def save(self, job: Job) -> None:
        await self._ensure_db()
        async with aiosqlite.connect(self.db_path) as db:
            await db.execute(
                """INSERT OR REPLACE INTO jobs
                   (id, tool_name, parameters, status, progress,
                    created_at, started_at, completed_at, result, error)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    job.id,
                    job.tool_name,
                    json.dumps(job.parameters),
                    job.status.value,
                    job.progress,
                    job.created_at.isoformat(),
                    job.started_at.isoformat() if job.started_at else None,
                    job.completed_at.isoformat() if job.completed_at else None,
                    json.dumps(job.result.to_dict()) if job.result else None,
                    job.error,
                ),
            )
            await db.commit()

    async def get(self, job_id: str) -> Job | None:
        await self._ensure_db()
        async with aiosqlite.connect(self.db_path) as db:
            db.row_factory = aiosqlite.Row
            cursor = await db.execute("SELECT * FROM jobs WHERE id = ? OR id LIKE ?", (job_id, f"{job_id}%"))
            row = await cursor.fetchone()
            if row is None:
                return None
            return self._row_to_job(row)

    async def list_all(self, status_filter: str | None = None) -> list[Job]:
        await self._ensure_db()
        async with aiosqlite.connect(self.db_path) as db:
            db.row_factory = aiosqlite.Row
            if status_filter:
                cursor = await db.execute(
                    "SELECT * FROM jobs WHERE status = ? ORDER BY created_at DESC", (status_filter,)
                )
            else:
                cursor = await db.execute("SELECT * FROM jobs ORDER BY created_at DESC")
            rows = await cursor.fetchall()
            return [self._row_to_job(row) for row in rows]

    def _row_to_job(self, row: Any) -> Job:
        result_data = json.loads(row["result"]) if row["result"] else None
        result = None
        if result_data:
            result = ToolResult(
                success=result_data["success"],
                data=result_data.get("data", {}),
                summary=result_data.get("summary", ""),
                warnings=result_data.get("warnings", []),
                error=result_data.get("error"),
            )

        return Job(
            id=row["id"],
            tool_name=row["tool_name"],
            parameters=json.loads(row["parameters"]),
            status=JobStatus(row["status"]),
            progress=row["progress"],
            created_at=datetime.fromisoformat(row["created_at"]),
            started_at=datetime.fromisoformat(row["started_at"]) if row["started_at"] else None,
            completed_at=datetime.fromisoformat(row["completed_at"]) if row["completed_at"] else None,
            result=result,
            error=row["error"],
        )
