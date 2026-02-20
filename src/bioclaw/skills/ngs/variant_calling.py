"""Variant calling tool using bcftools or GATK."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.ngs.variant_calling")


class VariantCallingTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="variant_calling",
            description=(
                "Call variants from an aligned BAM file using bcftools mpileup/call "
                "or GATK HaplotypeCaller. Produces a VCF file and reports summary "
                "statistics including total variants, SNPs, and indels."
            ),
            parameters=[
                ToolParameter(
                    name="bam_file",
                    type=ParameterType.STRING,
                    description="Path to the input sorted BAM file",
                ),
                ToolParameter(
                    name="reference",
                    type=ParameterType.STRING,
                    description="Path to the reference genome FASTA file",
                ),
                ToolParameter(
                    name="caller",
                    type=ParameterType.STRING,
                    description="Variant caller to use: 'bcftools' or 'gatk'",
                    required=False,
                    default="bcftools",
                    enum=["bcftools", "gatk"],
                ),
                ToolParameter(
                    name="regions",
                    type=ParameterType.STRING,
                    description=(
                        "Genomic regions to restrict calling (e.g., 'chr1:1000-2000' or "
                        "path to a BED file)"
                    ),
                    required=False,
                    default=None,
                ),
                ToolParameter(
                    name="output_vcf",
                    type=ParameterType.STRING,
                    description="Path for the output VCF file. Auto-generated if not specified.",
                    required=False,
                    default=None,
                ),
            ],
            category="ngs",
            is_long_running=True,
            required_binaries=[],  # checked dynamically based on caller
        )

    def check_available(self) -> bool:
        return shutil.which("bcftools") is not None or shutil.which("gatk") is not None

    async def execute(self, **kwargs) -> ToolResult:
        bam_file = kwargs["bam_file"]
        reference = kwargs["reference"]
        caller = kwargs.get("caller", "bcftools")
        regions = kwargs.get("regions")
        output_vcf = kwargs.get("output_vcf")

        # Validate input files
        for label, path_str in [("BAM file", bam_file), ("Reference", reference)]:
            if not Path(path_str).exists():
                return ToolResult(
                    success=False,
                    error=f"{label} not found: {path_str}",
                    summary=f"{label} does not exist: {path_str}",
                )

        # Check caller availability
        if caller == "bcftools" and shutil.which("bcftools") is None:
            return ToolResult(
                success=False,
                error="bcftools is not installed or not in PATH",
                summary="bcftools not available",
            )
        if caller == "gatk" and shutil.which("gatk") is None:
            return ToolResult(
                success=False,
                error="GATK is not installed or not in PATH",
                summary="GATK not available",
            )

        # Determine output VCF path
        if output_vcf is None:
            stem = Path(bam_file).stem.replace(".sorted", "").replace(".bam", "")
            output_vcf = str(Path(bam_file).parent / f"{stem}.vcf.gz")

        output_path = Path(output_vcf)
        os.makedirs(output_path.parent, exist_ok=True)

        try:
            if caller == "bcftools":
                return self._run_bcftools(bam_file, reference, regions, output_vcf)
            else:
                return self._run_gatk(bam_file, reference, regions, output_vcf)
        except subprocess.TimeoutExpired:
            return ToolResult(
                success=False,
                error="Variant calling timed out after 7200 seconds",
                summary="Variant calling timed out",
            )
        except Exception as e:
            return ToolResult(
                success=False,
                error=str(e),
                summary=f"Variant calling failed: {e}",
            )

    def _run_bcftools(
        self, bam_file: str, reference: str, regions: str | None, output_vcf: str
    ) -> ToolResult:
        # bcftools mpileup | bcftools call
        mpileup_cmd = [
            "bcftools", "mpileup",
            "-f", reference,
            "--max-depth", "250",
            "-Ou",
        ]
        if regions:
            if Path(regions).exists():
                mpileup_cmd.extend(["--regions-file", regions])
            else:
                mpileup_cmd.extend(["--regions", regions])
        mpileup_cmd.append(bam_file)

        call_cmd = ["bcftools", "call", "-mv", "-Oz", "-o", output_vcf]

        logger.info(f"Running: {' '.join(mpileup_cmd)} | {' '.join(call_cmd)}")

        mpileup_proc = subprocess.Popen(
            mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        call_proc = subprocess.Popen(
            call_cmd, stdin=mpileup_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        mpileup_proc.stdout.close()

        call_stdout, call_stderr = call_proc.communicate(timeout=7200)
        mpileup_stderr = mpileup_proc.stderr.read().decode()
        mpileup_proc.stderr.close()
        mpileup_proc.wait()

        if call_proc.returncode != 0:
            return ToolResult(
                success=False,
                error=f"bcftools call failed: {call_stderr.decode()}",
                summary="bcftools variant calling failed",
            )

        # Index VCF
        subprocess.run(
            ["bcftools", "index", output_vcf],
            capture_output=True, text=True, timeout=300,
        )

        return self._parse_vcf_stats(output_vcf, "bcftools")

    def _run_gatk(
        self, bam_file: str, reference: str, regions: str | None, output_vcf: str
    ) -> ToolResult:
        # Ensure reference dict and fai exist
        ref_path = Path(reference)
        dict_file = ref_path.with_suffix(".dict")
        fai_file = Path(reference + ".fai")

        if not fai_file.exists():
            result = subprocess.run(
                ["samtools", "faidx", reference],
                capture_output=True, text=True, timeout=600,
            )
            if result.returncode != 0:
                return ToolResult(
                    success=False,
                    error=f"samtools faidx failed: {result.stderr}",
                    summary="Reference indexing failed",
                )

        if not dict_file.exists() and shutil.which("samtools") is not None:
            result = subprocess.run(
                ["samtools", "dict", reference, "-o", str(dict_file)],
                capture_output=True, text=True, timeout=600,
            )

        # GATK HaplotypeCaller
        gatk_cmd = [
            "gatk", "HaplotypeCaller",
            "-R", reference,
            "-I", bam_file,
            "-O", output_vcf,
        ]
        if regions:
            if Path(regions).exists():
                gatk_cmd.extend(["-L", regions])
            else:
                gatk_cmd.extend(["-L", regions])

        logger.info(f"Running: {' '.join(gatk_cmd)}")

        result = subprocess.run(
            gatk_cmd,
            capture_output=True,
            text=True,
            timeout=7200,
        )

        if result.returncode != 0:
            return ToolResult(
                success=False,
                error=f"GATK HaplotypeCaller failed: {result.stderr[-2000:]}",
                summary="GATK variant calling failed",
            )

        return self._parse_vcf_stats(output_vcf, "gatk")

    def _parse_vcf_stats(self, vcf_path: str, caller: str) -> ToolResult:
        """Parse VCF to extract basic variant statistics."""
        total_variants = 0
        snps = 0
        indels = 0
        other = 0

        try:
            # Use bcftools stats if available for reliable parsing
            if shutil.which("bcftools"):
                stats_result = subprocess.run(
                    ["bcftools", "stats", vcf_path],
                    capture_output=True, text=True, timeout=300,
                )
                if stats_result.returncode == 0:
                    for line in stats_result.stdout.splitlines():
                        if line.startswith("SN\t") and "number of records:" in line:
                            total_variants = int(line.split("\t")[-1].strip())
                        elif line.startswith("SN\t") and "number of SNPs:" in line:
                            snps = int(line.split("\t")[-1].strip())
                        elif line.startswith("SN\t") and "number of indels:" in line:
                            indels = int(line.split("\t")[-1].strip())
                    other = total_variants - snps - indels
            else:
                # Fallback: count lines in VCF manually
                import gzip

                open_fn = gzip.open if vcf_path.endswith(".gz") else open
                with open_fn(vcf_path, "rt") as fh:
                    for line in fh:
                        if line.startswith("#"):
                            continue
                        total_variants += 1
                        parts = line.split("\t")
                        if len(parts) >= 5:
                            ref = parts[3]
                            alt = parts[4].split(",")[0]
                            if len(ref) == 1 and len(alt) == 1:
                                snps += 1
                            else:
                                indels += 1
                other = total_variants - snps - indels

        except Exception as e:
            logger.warning(f"VCF stats parsing failed: {e}")

        data = {
            "caller": caller,
            "output_vcf": vcf_path,
            "total_variants": total_variants,
            "snps": snps,
            "indels": indels,
            "other": other,
        }

        summary = (
            f"Variant calling with {caller} completed: "
            f"{total_variants:,} variants ({snps:,} SNPs, {indels:,} indels). "
            f"Output: {vcf_path}"
        )

        warnings = []
        if total_variants == 0:
            warnings.append("No variants called. Check input BAM quality and region coverage.")

        return ToolResult(
            success=True,
            data=data,
            summary=summary,
            files=[Path(vcf_path)],
            warnings=warnings,
        )
