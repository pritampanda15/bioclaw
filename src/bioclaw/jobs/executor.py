"""Job execution with process pool."""

from __future__ import annotations

import asyncio
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Callable

from bioclaw.utils.logging import get_logger

logger = get_logger("jobs.executor")

_pool: ProcessPoolExecutor | None = None


def get_process_pool(max_workers: int = 4) -> ProcessPoolExecutor:
    global _pool
    if _pool is None:
        _pool = ProcessPoolExecutor(max_workers=max_workers)
    return _pool


async def run_in_process(func: Callable, *args: Any, **kwargs: Any) -> Any:
    """Run a CPU-bound function in the process pool."""
    loop = asyncio.get_event_loop()
    pool = get_process_pool()
    return await loop.run_in_executor(pool, func, *args)


async def run_subprocess(
    cmd: list[str],
    cwd: str | None = None,
    timeout: int = 3600,
) -> tuple[int, str, str]:
    """Run an external command as async subprocess.

    Returns (returncode, stdout, stderr).
    """
    logger.debug(f"Running subprocess: {' '.join(cmd)}")
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=cwd,
    )
    try:
        stdout, stderr = await asyncio.wait_for(proc.communicate(), timeout=timeout)
        return (
            proc.returncode or 0,
            stdout.decode("utf-8", errors="replace"),
            stderr.decode("utf-8", errors="replace"),
        )
    except asyncio.TimeoutError:
        proc.kill()
        await proc.wait()
        raise TimeoutError(f"Subprocess timed out after {timeout}s: {' '.join(cmd)}")
