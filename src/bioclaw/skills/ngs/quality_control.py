"""FastQC quality control wrapper tool."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.ngs.quality_control")


class QualityControlTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="quality_control",
            description=(
                "Run FastQC quality control analysis on a FASTQ file. "
                "Returns per-module pass/warn/fail status, basic statistics "
                "(total sequences, sequence length, GC content), and the path "
                "to the full HTML report."
            ),
            parameters=[
                ToolParameter(
                    name="input_file",
                    type=ParameterType.STRING,
                    description="Path to the input FASTQ file (.fastq or .fastq.gz)",
                ),
                ToolParameter(
                    name="output_dir",
                    type=ParameterType.STRING,
                    description="Directory for FastQC output. Defaults to a temp directory.",
                    required=False,
                    default=None,
                ),
            ],
            category="ngs",
            is_long_running=True,
            required_binaries=["fastqc"],
        )

    def check_available(self) -> bool:
        return shutil.which("fastqc") is not None

    async def execute(self, **kwargs) -> ToolResult:
        input_file = kwargs["input_file"]
        output_dir = kwargs.get("output_dir")

        # Validate input file
        input_path = Path(input_file)
        if not input_path.exists():
            return ToolResult(
                success=False,
                error=f"Input file not found: {input_file}",
                summary=f"FASTQ file does not exist: {input_file}",
            )

        # Create output directory
        temp_dir = None
        if output_dir is None:
            temp_dir = tempfile.mkdtemp(prefix="bioclaw_fastqc_")
            output_dir = temp_dir
        else:
            os.makedirs(output_dir, exist_ok=True)

        try:
            # Run FastQC
            cmd = ["fastqc", str(input_path), "-o", output_dir, "--extract", "--quiet"]
            logger.info(f"Running FastQC: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600,
            )

            if result.returncode != 0:
                return ToolResult(
                    success=False,
                    error=f"FastQC failed (exit code {result.returncode}): {result.stderr}",
                    summary="FastQC execution failed",
                )

            # Locate the extracted output directory
            stem = input_path.name.replace(".fastq.gz", "").replace(".fastq", "")
            fastqc_dir = Path(output_dir) / f"{stem}_fastqc"
            summary_file = fastqc_dir / "summary.txt"
            fastqc_data_file = fastqc_dir / "fastqc_data.txt"
            html_report = Path(output_dir) / f"{stem}_fastqc.html"

            # Parse summary.txt for module pass/warn/fail
            modules: dict[str, str] = {}
            if summary_file.exists():
                for line in summary_file.read_text().splitlines():
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        status = parts[0]  # PASS, WARN, or FAIL
                        module_name = parts[1]
                        modules[module_name] = status

            # Parse basic statistics from fastqc_data.txt
            basic_stats: dict[str, str] = {}
            if fastqc_data_file.exists():
                in_basic_stats = False
                for line in fastqc_data_file.read_text().splitlines():
                    if line.startswith(">>Basic Statistics"):
                        in_basic_stats = True
                        continue
                    if in_basic_stats and line.startswith(">>END_MODULE"):
                        break
                    if in_basic_stats and not line.startswith("#"):
                        parts = line.split("\t")
                        if len(parts) >= 2:
                            basic_stats[parts[0]] = parts[1]

            pass_count = sum(1 for v in modules.values() if v == "PASS")
            warn_count = sum(1 for v in modules.values() if v == "WARN")
            fail_count = sum(1 for v in modules.values() if v == "FAIL")

            output_files = []
            if html_report.exists():
                output_files.append(html_report)

            data = {
                "modules": modules,
                "basic_statistics": basic_stats,
                "pass_count": pass_count,
                "warn_count": warn_count,
                "fail_count": fail_count,
                "html_report": str(html_report) if html_report.exists() else None,
                "output_dir": output_dir,
            }

            warnings = []
            if warn_count > 0:
                warned = [m for m, s in modules.items() if s == "WARN"]
                warnings.append(f"Warnings in: {', '.join(warned)}")
            if fail_count > 0:
                failed = [m for m, s in modules.items() if s == "FAIL"]
                warnings.append(f"Failures in: {', '.join(failed)}")

            total_seqs = basic_stats.get("Total Sequences", "N/A")
            seq_length = basic_stats.get("Sequence length", "N/A")
            gc_content = basic_stats.get("%GC", "N/A")

            summary = (
                f"FastQC completed for {input_path.name}: "
                f"{total_seqs} sequences, length {seq_length}, GC {gc_content}%. "
                f"Modules: {pass_count} passed, {warn_count} warnings, {fail_count} failed."
            )

            return ToolResult(
                success=True,
                data=data,
                summary=summary,
                files=output_files,
                warnings=warnings,
            )

        except subprocess.TimeoutExpired:
            return ToolResult(
                success=False,
                error="FastQC timed out after 3600 seconds",
                summary="FastQC execution timed out",
            )
        except Exception as e:
            return ToolResult(
                success=False,
                error=str(e),
                summary=f"FastQC failed: {e}",
            )
