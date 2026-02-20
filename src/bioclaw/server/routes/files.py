"""File upload/download routes."""

from __future__ import annotations

import shutil
from pathlib import Path

from fastapi import APIRouter, HTTPException, UploadFile
from fastapi.responses import FileResponse

from bioclaw.config.settings import get_settings

router = APIRouter(tags=["files"])


@router.post("/files/upload")
async def upload_file(file: UploadFile) -> dict:
    """Upload a file to the BioClaw data directory."""
    settings = get_settings()
    uploads_dir = settings.data_dir / "uploads"
    uploads_dir.mkdir(parents=True, exist_ok=True)

    if not file.filename:
        raise HTTPException(status_code=400, detail="No filename provided")

    # Sanitize filename
    safe_name = Path(file.filename).name
    if not safe_name or safe_name.startswith("."):
        raise HTTPException(status_code=400, detail="Invalid filename")

    dest = uploads_dir / safe_name
    with open(dest, "wb") as f:
        shutil.copyfileobj(file.file, f)

    return {"filename": safe_name, "path": str(dest), "size": dest.stat().st_size}


@router.get("/files/{filename}")
async def download_file(filename: str) -> FileResponse:
    """Download a file from results or uploads."""
    settings = get_settings()

    # Search in uploads and results
    for subdir in ["uploads", "results"]:
        search_dir = settings.data_dir / subdir
        if not search_dir.exists():
            continue
        # Direct match
        candidate = search_dir / filename
        if candidate.exists() and candidate.is_file():
            return FileResponse(candidate, filename=filename)
        # Recursive search
        for f in search_dir.rglob(filename):
            if f.is_file():
                return FileResponse(f, filename=filename)

    raise HTTPException(status_code=404, detail=f"File not found: {filename}")


@router.get("/files")
async def list_files() -> list[dict]:
    """List files in the data directory."""
    settings = get_settings()
    files = []
    for subdir in ["uploads", "results"]:
        search_dir = settings.data_dir / subdir
        if not search_dir.exists():
            continue
        for f in search_dir.rglob("*"):
            if f.is_file():
                files.append({
                    "name": f.name,
                    "path": str(f.relative_to(settings.data_dir)),
                    "size": f.stat().st_size,
                })
    return files
