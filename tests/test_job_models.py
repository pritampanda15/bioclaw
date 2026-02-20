"""Tests for job models."""

from bioclaw.jobs.models import Job, JobStatus


def test_job_creation():
    job = Job(tool_name="test_tool", parameters={"key": "value"})
    assert job.status == JobStatus.PENDING
    assert job.tool_name == "test_tool"
    assert len(job.id) == 32  # uuid hex


def test_job_to_dict():
    job = Job(tool_name="test_tool", parameters={"smiles": "CCO"})
    d = job.to_dict()
    assert d["tool_name"] == "test_tool"
    assert d["status"] == "pending"
    assert d["parameters"]["smiles"] == "CCO"
    assert "id" in d
    assert "created_at" in d
