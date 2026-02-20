"""Job management API routes."""

from __future__ import annotations

from fastapi import APIRouter, HTTPException

from bioclaw.jobs.manager import get_job_manager

router = APIRouter(tags=["jobs"])


@router.get("/jobs")
async def list_jobs(status: str | None = None) -> list[dict]:
    """List all jobs, optionally filtered by status."""
    manager = get_job_manager()
    jobs = await manager.list_jobs(status_filter=status)
    return [job.to_dict() for job in jobs]


@router.get("/jobs/{job_id}")
async def get_job(job_id: str) -> dict:
    """Get a specific job by ID."""
    manager = get_job_manager()
    job = await manager.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")
    return job.to_dict()


@router.post("/jobs/{job_id}/cancel")
async def cancel_job(job_id: str) -> dict:
    """Cancel a running job."""
    manager = get_job_manager()
    success = await manager.cancel_job(job_id)
    if not success:
        raise HTTPException(status_code=400, detail="Failed to cancel job")
    return {"status": "cancelled", "job_id": job_id}
