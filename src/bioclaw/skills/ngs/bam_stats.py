"""BAM file statistics tool using pysam or samtools."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.ngs.bam_stats")


def _pysam_available() -> bool:
    try:
        import pysam  # noqa: F401

        return True
    except ImportError:
        return False


class BamStatsTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="bam_stats",
            description=(
                "Compute statistics for a BAM file: total reads, mapped reads, "
                "mapping rate, mean coverage, and insert size statistics "
                "(mean, std dev). Uses pysam when available, falling back to "
                "samtools if pysam is not installed."
            ),
            parameters=[
                ToolParameter(
                    name="bam_file",
                    type=ParameterType.STRING,
                    description="Path to the input BAM file",
                ),
            ],
            category="ngs",
            required_binaries=[],  # samtools only needed as fallback
        )

    def check_available(self) -> bool:
        if _pysam_available():
            return True
        if shutil.which("samtools") is not None:
            return True
        return False

    async def execute(self, **kwargs) -> ToolResult:
        bam_file = kwargs["bam_file"]

        bam_path = Path(bam_file)
        if not bam_path.exists():
            return ToolResult(
                success=False,
                error=f"BAM file not found: {bam_file}",
                summary=f"BAM file does not exist: {bam_file}",
            )

        if _pysam_available():
            return self._stats_with_pysam(bam_path)
        elif shutil.which("samtools") is not None:
            return self._stats_with_samtools(bam_path)
        else:
            return ToolResult(
                success=False,
                error="Neither pysam nor samtools is available",
                summary="No BAM processing backend available. Install pysam or samtools.",
            )

    def _stats_with_pysam(self, bam_path: Path) -> ToolResult:
        try:
            import pysam
        except ImportError:
            return ToolResult(
                success=False,
                error="pysam import failed unexpectedly",
                summary="pysam could not be imported",
            )

        try:
            bam = pysam.AlignmentFile(str(bam_path), "rb")
        except Exception as e:
            return ToolResult(
                success=False,
                error=f"Failed to open BAM file: {e}",
                summary=f"Could not open BAM: {e}",
            )

        try:
            total_reads = 0
            mapped_reads = 0
            insert_sizes: list[int] = []

            # Coverage tracking
            coverage_bases = 0
            genome_length = 0

            # Get reference lengths for coverage calculation
            for length in bam.lengths:
                genome_length += length

            for read in bam.fetch(until_eof=True):
                total_reads += 1
                if not read.is_unmapped:
                    mapped_reads += 1
                    if read.query_alignment_length:
                        coverage_bases += read.query_alignment_length
                    # Insert size for properly paired reads
                    if (
                        read.is_proper_pair
                        and not read.is_secondary
                        and not read.is_supplementary
                        and read.template_length > 0
                    ):
                        insert_sizes.append(read.template_length)

            bam.close()

            mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0.0
            mean_coverage = (coverage_bases / genome_length) if genome_length > 0 else 0.0

            # Insert size statistics
            insert_size_stats: dict[str, float | None] = {}
            if insert_sizes:
                mean_insert = sum(insert_sizes) / len(insert_sizes)
                variance = sum((x - mean_insert) ** 2 for x in insert_sizes) / len(insert_sizes)
                std_insert = variance ** 0.5
                insert_size_stats = {
                    "mean": round(mean_insert, 2),
                    "std_dev": round(std_insert, 2),
                    "min": min(insert_sizes),
                    "max": max(insert_sizes),
                    "count": len(insert_sizes),
                }
            else:
                insert_size_stats = {
                    "mean": None,
                    "std_dev": None,
                    "min": None,
                    "max": None,
                    "count": 0,
                }

            data = {
                "bam_file": str(bam_path),
                "backend": "pysam",
                "total_reads": total_reads,
                "mapped_reads": mapped_reads,
                "unmapped_reads": total_reads - mapped_reads,
                "mapping_rate": round(mapping_rate, 2),
                "mean_coverage": round(mean_coverage, 2),
                "genome_length": genome_length,
                "insert_size": insert_size_stats,
            }

            summary = (
                f"BAM stats for {bam_path.name}: "
                f"{total_reads:,} total reads, {mapped_reads:,} mapped "
                f"({mapping_rate:.1f}%), mean coverage {mean_coverage:.1f}x"
            )
            if insert_size_stats.get("mean") is not None:
                summary += f", mean insert size {insert_size_stats['mean']:.0f} bp"

            return ToolResult(success=True, data=data, summary=summary)

        except Exception as e:
            return ToolResult(
                success=False,
                error=str(e),
                summary=f"Error computing BAM stats with pysam: {e}",
            )

    def _stats_with_samtools(self, bam_path: Path) -> ToolResult:
        try:
            # Run samtools flagstat
            flagstat_result = subprocess.run(
                ["samtools", "flagstat", str(bam_path)],
                capture_output=True,
                text=True,
                timeout=600,
            )
            if flagstat_result.returncode != 0:
                return ToolResult(
                    success=False,
                    error=f"samtools flagstat failed: {flagstat_result.stderr}",
                    summary="samtools flagstat failed",
                )

            # Parse flagstat output
            total_reads = 0
            mapped_reads = 0
            for line in flagstat_result.stdout.splitlines():
                if "in total" in line:
                    total_reads = int(line.split()[0])
                elif "mapped (" in line and "primary" not in line:
                    mapped_reads = int(line.split()[0])

            mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0.0

            # Run samtools stats for more detailed metrics
            stats_result = subprocess.run(
                ["samtools", "stats", str(bam_path)],
                capture_output=True,
                text=True,
                timeout=600,
            )

            mean_coverage = 0.0
            insert_size_mean = None
            insert_size_std = None
            bases_mapped = 0
            genome_length = 0

            if stats_result.returncode == 0:
                for line in stats_result.stdout.splitlines():
                    if not line.startswith("SN\t"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue
                    key = parts[1].rstrip(":")
                    value = parts[2].strip()
                    if key == "bases mapped (cigar)":
                        bases_mapped = int(value)
                    elif key == "total length":
                        genome_length = int(value) if int(value) > 0 else genome_length
                    elif key == "insert size average":
                        insert_size_mean = float(value)
                    elif key == "insert size standard deviation":
                        insert_size_std = float(value)

                if genome_length > 0:
                    mean_coverage = bases_mapped / genome_length

            insert_size_stats: dict[str, float | None] = {
                "mean": round(insert_size_mean, 2) if insert_size_mean is not None else None,
                "std_dev": round(insert_size_std, 2) if insert_size_std is not None else None,
            }

            data = {
                "bam_file": str(bam_path),
                "backend": "samtools",
                "total_reads": total_reads,
                "mapped_reads": mapped_reads,
                "unmapped_reads": total_reads - mapped_reads,
                "mapping_rate": round(mapping_rate, 2),
                "mean_coverage": round(mean_coverage, 2),
                "genome_length": genome_length,
                "insert_size": insert_size_stats,
            }

            summary = (
                f"BAM stats for {bam_path.name} (samtools): "
                f"{total_reads:,} total reads, {mapped_reads:,} mapped "
                f"({mapping_rate:.1f}%), mean coverage {mean_coverage:.1f}x"
            )
            if insert_size_mean is not None:
                summary += f", mean insert size {insert_size_mean:.0f} bp"

            return ToolResult(success=True, data=data, summary=summary)

        except subprocess.TimeoutExpired:
            return ToolResult(
                success=False,
                error="samtools timed out after 600 seconds",
                summary="samtools execution timed out",
            )
        except Exception as e:
            return ToolResult(
                success=False,
                error=str(e),
                summary=f"Error computing BAM stats with samtools: {e}",
            )
