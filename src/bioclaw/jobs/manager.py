"""Job lifecycle manager."""

from __future__ import annotations

import asyncio
from datetime import datetime
from typing import Any

from bioclaw.config.settings import get_settings
from bioclaw.jobs.models import Job, JobStatus
from bioclaw.jobs.store import JobStore
from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("jobs.manager")


class JobManager:
    """Manages async job lifecycle with SQLite persistence."""

    def __init__(self, store: JobStore) -> None:
        self.store = store
        self._running_tasks: dict[str, asyncio.Task] = {}

    async def submit(self, tool: BaseTool, parameters: dict[str, Any]) -> Job:
        """Submit a new job for async execution."""
        job = Job(tool_name=tool.metadata.name, parameters=parameters)
        await self.store.save(job)

        task = asyncio.create_task(self._execute(job, tool))
        self._running_tasks[job.id] = task
        logger.info(f"Job {job.id[:8]} submitted: {tool.metadata.name}")
        return job

    async def _execute(self, job: Job, tool: BaseTool) -> None:
        """Execute a job and update its status."""
        job.status = JobStatus.RUNNING
        job.started_at = datetime.now()
        await self.store.save(job)

        try:
            result = await tool.safe_execute(**job.parameters)
            job.status = JobStatus.COMPLETED if result.success else JobStatus.FAILED
            job.result = result
            if not result.success:
                job.error = result.error
        except asyncio.CancelledError:
            job.status = JobStatus.CANCELLED
            job.result = ToolResult(success=False, error="Job cancelled", summary="Job was cancelled")
        except Exception as e:
            logger.exception(f"Job {job.id[:8]} failed")
            job.status = JobStatus.FAILED
            job.error = str(e)
            job.result = ToolResult(success=False, error=str(e), summary=f"Job failed: {e}")
        finally:
            job.completed_at = datetime.now()
            await self.store.save(job)
            self._running_tasks.pop(job.id, None)
            logger.info(f"Job {job.id[:8]} {job.status.value}")

    async def get_job(self, job_id: str) -> Job | None:
        return await self.store.get(job_id)

    async def list_jobs(self, status_filter: str | None = None) -> list[Job]:
        return await self.store.list_all(status_filter)

    async def cancel_job(self, job_id: str) -> bool:
        job = await self.store.get(job_id)
        if job is None:
            return False
        task = self._running_tasks.get(job.id)
        if task and not task.done():
            task.cancel()
            return True
        return False

    async def get_job_status_for_agent(self, job_id: str) -> ToolResult:
        """Get job status formatted for the agent to consume."""
        job = await self.get_job(job_id)
        if job is None:
            return ToolResult(success=False, error=f"Job {job_id} not found")
        data = job.to_dict()
        summary = f"Job {job.id[:8]}: {job.status.value}"
        if job.result:
            summary += f" — {job.result.summary}"
        return ToolResult(success=True, data=data, summary=summary)


_manager: JobManager | None = None


def get_job_manager() -> JobManager:
    global _manager
    if _manager is None:
        settings = get_settings()
        store = JobStore(settings.data_dir / "jobs.db")
        _manager = JobManager(store)
    return _manager
