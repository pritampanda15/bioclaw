"""SQLite database setup and utilities."""

from __future__ import annotations

from pathlib import Path

import aiosqlite

from bioclaw.config.settings import get_settings
from bioclaw.utils.logging import get_logger

logger = get_logger("storage.database")


async def get_db_connection(db_name: str = "bioclaw.db") -> aiosqlite.Connection:
    """Get a connection to the SQLite database."""
    settings = get_settings()
    db_path = settings.data_dir / db_name
    db_path.parent.mkdir(parents=True, exist_ok=True)
    return await aiosqlite.connect(db_path)


async def init_database() -> None:
    """Initialize all database tables."""
    settings = get_settings()
    db_path = settings.data_dir / "bioclaw.db"
    db_path.parent.mkdir(parents=True, exist_ok=True)

    async with aiosqlite.connect(db_path) as db:
        await db.execute("""
            CREATE TABLE IF NOT EXISTS chat_sessions (
                id TEXT PRIMARY KEY,
                title TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
        """)
        await db.execute("""
            CREATE TABLE IF NOT EXISTS chat_messages (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                session_id TEXT NOT NULL,
                role TEXT NOT NULL,
                content TEXT NOT NULL,
                created_at TEXT NOT NULL,
                FOREIGN KEY (session_id) REFERENCES chat_sessions(id)
            )
        """)
        await db.execute("""
            CREATE TABLE IF NOT EXISTS results (
                id TEXT PRIMARY KEY,
                job_id TEXT,
                tool_name TEXT NOT NULL,
                file_path TEXT,
                data TEXT,
                created_at TEXT NOT NULL
            )
        """)
        await db.commit()
    logger.info(f"Database initialized at {db_path}")
