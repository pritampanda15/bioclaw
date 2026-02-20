"""Read alignment tool using BWA or HISAT2."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.ngs.read_alignment")


class ReadAlignmentTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="read_alignment",
            description=(
                "Align sequencing reads to a reference genome using BWA-MEM or HISAT2. "
                "Supports single-end and paired-end reads. Output is a coordinate-sorted "
                "BAM file. Returns alignment statistics including total reads, mapping rate, "
                "and properly paired percentage."
            ),
            parameters=[
                ToolParameter(
                    name="reads_1",
                    type=ParameterType.STRING,
                    description="Path to the first FASTQ file (R1)",
                ),
                ToolParameter(
                    name="reads_2",
                    type=ParameterType.STRING,
                    description="Path to the second FASTQ file (R2) for paired-end reads",
                    required=False,
                    default=None,
                ),
                ToolParameter(
                    name="reference",
                    type=ParameterType.STRING,
                    description="Path to the reference genome FASTA file",
                ),
                ToolParameter(
                    name="aligner",
                    type=ParameterType.STRING,
                    description="Aligner to use: 'bwa' or 'hisat2'",
                    required=False,
                    default="bwa",
                    enum=["bwa", "hisat2"],
                ),
                ToolParameter(
                    name="threads",
                    type=ParameterType.INTEGER,
                    description="Number of threads to use",
                    required=False,
                    default=4,
                ),
                ToolParameter(
                    name="output_bam",
                    type=ParameterType.STRING,
                    description="Path for the output sorted BAM file. Auto-generated if not specified.",
                    required=False,
                    default=None,
                ),
            ],
            category="ngs",
            is_long_running=True,
            required_binaries=["samtools"],  # aligner checked dynamically
        )

    def check_available(self) -> bool:
        if shutil.which("samtools") is None:
            return False
        # At least one aligner must be present
        return shutil.which("bwa") is not None or shutil.which("hisat2") is not None

    async def execute(self, **kwargs) -> ToolResult:
        reads_1 = kwargs["reads_1"]
        reads_2 = kwargs.get("reads_2")
        reference = kwargs["reference"]
        aligner = kwargs.get("aligner", "bwa")
        threads = kwargs.get("threads", 4)
        output_bam = kwargs.get("output_bam")

        # Validate input files
        for label, path_str in [("Reads R1", reads_1), ("Reference", reference)]:
            if not Path(path_str).exists():
                return ToolResult(
                    success=False,
                    error=f"{label} file not found: {path_str}",
                    summary=f"{label} file does not exist: {path_str}",
                )
        if reads_2 and not Path(reads_2).exists():
            return ToolResult(
                success=False,
                error=f"Reads R2 file not found: {reads_2}",
                summary=f"Reads R2 file does not exist: {reads_2}",
            )

        # Check selected aligner is available
        if shutil.which(aligner) is None:
            return ToolResult(
                success=False,
                error=f"Aligner not found: {aligner}",
                summary=f"{aligner} is not installed or not in PATH",
            )

        # Determine output BAM path
        if output_bam is None:
            stem = Path(reads_1).stem.replace(".fastq", "").replace(".fq", "")
            output_bam = str(Path(reads_1).parent / f"{stem}.sorted.bam")

        output_path = Path(output_bam)
        os.makedirs(output_path.parent, exist_ok=True)

        try:
            if aligner == "bwa":
                return await self._run_bwa(reads_1, reads_2, reference, threads, output_bam)
            else:
                return await self._run_hisat2(reads_1, reads_2, reference, threads, output_bam)
        except subprocess.TimeoutExpired:
            return ToolResult(
                success=False,
                error="Alignment timed out after 7200 seconds",
                summary="Alignment execution timed out",
            )
        except Exception as e:
            return ToolResult(
                success=False,
                error=str(e),
                summary=f"Alignment failed: {e}",
            )

    async def _run_bwa(
        self, reads_1: str, reads_2: str | None, reference: str, threads: int, output_bam: str
    ) -> ToolResult:
        # Check if BWA index exists; if not, build it
        index_extensions = [".bwt", ".pac", ".ann", ".amb", ".sa"]
        if not all(Path(reference + ext).exists() for ext in index_extensions):
            logger.info("BWA index not found, building index...")
            idx_result = subprocess.run(
                ["bwa", "index", reference],
                capture_output=True,
                text=True,
                timeout=7200,
            )
            if idx_result.returncode != 0:
                return ToolResult(
                    success=False,
                    error=f"BWA indexing failed: {idx_result.stderr}",
                    summary="BWA index build failed",
                )

        # Build BWA MEM command
        bwa_cmd = ["bwa", "mem", "-t", str(threads), reference, reads_1]
        if reads_2:
            bwa_cmd.append(reads_2)

        # Pipe: bwa mem | samtools sort
        sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", output_bam, "-"]

        logger.info(f"Running: {' '.join(bwa_cmd)} | {' '.join(sort_cmd)}")

        bwa_proc = subprocess.Popen(
            bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        sort_proc = subprocess.Popen(
            sort_cmd, stdin=bwa_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        bwa_proc.stdout.close()  # Allow bwa to receive SIGPIPE if sort exits

        sort_stdout, sort_stderr = sort_proc.communicate(timeout=7200)
        bwa_stderr = bwa_proc.stderr.read().decode()
        bwa_proc.stderr.close()
        bwa_proc.wait()

        if sort_proc.returncode != 0:
            return ToolResult(
                success=False,
                error=f"samtools sort failed: {sort_stderr.decode()}",
                summary="BAM sorting failed",
            )

        # Index the BAM
        idx_result = subprocess.run(
            ["samtools", "index", output_bam],
            capture_output=True,
            text=True,
            timeout=600,
        )

        return self._collect_stats(output_bam, "bwa", reads_2 is not None, bwa_stderr)

    async def _run_hisat2(
        self, reads_1: str, reads_2: str | None, reference: str, threads: int, output_bam: str
    ) -> ToolResult:
        # Check if HISAT2 index exists
        index_base = reference.rsplit(".", 1)[0] if "." in reference else reference
        hisat2_index_files = [f"{index_base}.1.ht2", f"{index_base}.1.ht2l"]
        index_exists = any(Path(f).exists() for f in hisat2_index_files)

        if not index_exists:
            logger.info("HISAT2 index not found, building index...")
            idx_result = subprocess.run(
                ["hisat2-build", reference, index_base],
                capture_output=True,
                text=True,
                timeout=7200,
            )
            if idx_result.returncode != 0:
                return ToolResult(
                    success=False,
                    error=f"HISAT2 index build failed: {idx_result.stderr}",
                    summary="HISAT2 index build failed",
                )

        # Build HISAT2 command
        hisat2_cmd = ["hisat2", "-x", index_base, "-p", str(threads)]
        if reads_2:
            hisat2_cmd.extend(["-1", reads_1, "-2", reads_2])
        else:
            hisat2_cmd.extend(["-U", reads_1])

        sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", output_bam, "-"]

        logger.info(f"Running: {' '.join(hisat2_cmd)} | {' '.join(sort_cmd)}")

        hisat2_proc = subprocess.Popen(
            hisat2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        sort_proc = subprocess.Popen(
            sort_cmd, stdin=hisat2_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        hisat2_proc.stdout.close()

        sort_stdout, sort_stderr = sort_proc.communicate(timeout=7200)
        hisat2_stderr = hisat2_proc.stderr.read().decode()
        hisat2_proc.stderr.close()
        hisat2_proc.wait()

        if sort_proc.returncode != 0:
            return ToolResult(
                success=False,
                error=f"samtools sort failed: {sort_stderr.decode()}",
                summary="BAM sorting failed",
            )

        # Index the BAM
        subprocess.run(
            ["samtools", "index", output_bam],
            capture_output=True,
            text=True,
            timeout=600,
        )

        return self._collect_stats(output_bam, "hisat2", reads_2 is not None, hisat2_stderr)

    def _collect_stats(
        self, output_bam: str, aligner: str, paired: bool, aligner_stderr: str
    ) -> ToolResult:
        """Run samtools flagstat and return alignment statistics."""
        flagstat = subprocess.run(
            ["samtools", "flagstat", output_bam],
            capture_output=True,
            text=True,
            timeout=300,
        )

        total_reads = 0
        mapped_reads = 0
        properly_paired = 0

        if flagstat.returncode == 0:
            for line in flagstat.stdout.splitlines():
                if "in total" in line:
                    total_reads = int(line.split()[0])
                elif "mapped (" in line and "primary" not in line:
                    mapped_reads = int(line.split()[0])
                elif "properly paired" in line:
                    properly_paired = int(line.split()[0])

        mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0.0
        paired_rate = (properly_paired / total_reads * 100) if total_reads > 0 else 0.0

        data = {
            "aligner": aligner,
            "output_bam": output_bam,
            "paired_end": paired,
            "total_reads": total_reads,
            "mapped_reads": mapped_reads,
            "mapping_rate": round(mapping_rate, 2),
            "properly_paired": properly_paired,
            "properly_paired_rate": round(paired_rate, 2),
            "aligner_log": aligner_stderr[-2000:] if aligner_stderr else "",
        }

        warnings = []
        if mapping_rate < 70:
            warnings.append(f"Low mapping rate: {mapping_rate:.1f}%")

        summary = (
            f"Alignment completed with {aligner}: {total_reads:,} reads, "
            f"{mapped_reads:,} mapped ({mapping_rate:.1f}%)"
        )
        if paired:
            summary += f", {properly_paired:,} properly paired ({paired_rate:.1f}%)"
        summary += f". Output: {output_bam}"

        return ToolResult(
            success=True,
            data=data,
            summary=summary,
            files=[Path(output_bam)],
            warnings=warnings,
        )
