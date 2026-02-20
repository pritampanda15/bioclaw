"""RNA-seq meta-pipeline: QC -> Alignment -> Feature Counting."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.ngs.rnaseq_pipeline")


class RnaseqPipelineTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="rnaseq_pipeline",
            description=(
                "Run a complete RNA-seq analysis pipeline: "
                "FastQC quality control -> HISAT2 alignment -> featureCounts gene quantification. "
                "Returns summary results from each step including QC metrics, alignment rates, "
                "and gene count statistics."
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
                    name="annotation_gtf",
                    type=ParameterType.STRING,
                    description="Path to the gene annotation GTF file",
                ),
                ToolParameter(
                    name="output_dir",
                    type=ParameterType.STRING,
                    description="Directory for all pipeline outputs",
                ),
                ToolParameter(
                    name="threads",
                    type=ParameterType.INTEGER,
                    description="Number of threads to use",
                    required=False,
                    default=4,
                ),
            ],
            category="ngs",
            is_long_running=True,
            required_binaries=["fastqc", "hisat2", "samtools", "featureCounts"],
        )

    def check_available(self) -> bool:
        for binary in self.metadata.required_binaries:
            if shutil.which(binary) is None:
                logger.debug(f"Binary not found: {binary}")
                return False
        return True

    async def execute(self, **kwargs) -> ToolResult:
        reads_1 = kwargs["reads_1"]
        reads_2 = kwargs.get("reads_2")
        reference = kwargs["reference"]
        annotation_gtf = kwargs["annotation_gtf"]
        output_dir = kwargs["output_dir"]
        threads = kwargs.get("threads", 4)

        # Validate input files
        required_files = [
            ("Reads R1", reads_1),
            ("Reference", reference),
            ("Annotation GTF", annotation_gtf),
        ]
        if reads_2:
            required_files.append(("Reads R2", reads_2))

        for label, path_str in required_files:
            if not Path(path_str).exists():
                return ToolResult(
                    success=False,
                    error=f"{label} not found: {path_str}",
                    summary=f"{label} does not exist: {path_str}",
                )

        # Create output directory structure
        os.makedirs(output_dir, exist_ok=True)
        qc_dir = os.path.join(output_dir, "fastqc")
        align_dir = os.path.join(output_dir, "alignment")
        counts_dir = os.path.join(output_dir, "counts")
        for d in [qc_dir, align_dir, counts_dir]:
            os.makedirs(d, exist_ok=True)

        step_results: dict[str, dict] = {}
        warnings: list[str] = []
        output_files: list[Path] = []

        # ---- Step 1: FastQC ----
        logger.info("Step 1/3: Running FastQC")
        try:
            qc_result = self._run_fastqc(reads_1, reads_2, qc_dir)
            step_results["fastqc"] = qc_result
            if qc_result.get("warnings"):
                warnings.extend(
                    [f"FastQC: {w}" for w in qc_result["warnings"]]
                )
        except Exception as e:
            step_results["fastqc"] = {"status": "failed", "error": str(e)}
            warnings.append(f"FastQC failed: {e}")

        # ---- Step 2: HISAT2 Alignment ----
        logger.info("Step 2/3: Running HISAT2 alignment")
        stem = Path(reads_1).stem.replace(".fastq", "").replace(".fq", "")
        output_bam = os.path.join(align_dir, f"{stem}.sorted.bam")

        try:
            align_result = self._run_hisat2(
                reads_1, reads_2, reference, threads, output_bam
            )
            step_results["alignment"] = align_result
            if align_result.get("status") == "success":
                output_files.append(Path(output_bam))
                if align_result.get("mapping_rate", 100) < 70:
                    warnings.append(
                        f"Low mapping rate: {align_result['mapping_rate']:.1f}%"
                    )
            else:
                return ToolResult(
                    success=False,
                    data={"steps": step_results},
                    error=f"Alignment failed: {align_result.get('error', 'unknown')}",
                    summary="Pipeline stopped: alignment step failed",
                    warnings=warnings,
                )
        except Exception as e:
            return ToolResult(
                success=False,
                data={"steps": step_results},
                error=f"Alignment failed: {e}",
                summary=f"Pipeline stopped at alignment: {e}",
                warnings=warnings,
            )

        # ---- Step 3: featureCounts ----
        logger.info("Step 3/3: Running featureCounts")
        counts_file = os.path.join(counts_dir, f"{stem}.counts.txt")

        try:
            counts_result = self._run_featurecounts(
                output_bam, annotation_gtf, counts_file, threads, reads_2 is not None
            )
            step_results["featurecounts"] = counts_result
            if counts_result.get("status") == "success":
                output_files.append(Path(counts_file))
            else:
                warnings.append(
                    f"featureCounts issues: {counts_result.get('error', 'unknown')}"
                )
        except Exception as e:
            step_results["featurecounts"] = {"status": "failed", "error": str(e)}
            warnings.append(f"featureCounts failed: {e}")

        # Build summary
        qc_status = step_results.get("fastqc", {}).get("status", "unknown")
        align_info = step_results.get("alignment", {})
        counts_info = step_results.get("featurecounts", {})

        summary_parts = [f"RNA-seq pipeline completed for {Path(reads_1).name}."]
        summary_parts.append(f"FastQC: {qc_status}.")
        if align_info.get("status") == "success":
            summary_parts.append(
                f"Alignment: {align_info.get('mapped_reads', 0):,} mapped "
                f"({align_info.get('mapping_rate', 0):.1f}%)."
            )
        if counts_info.get("status") == "success":
            summary_parts.append(
                f"featureCounts: {counts_info.get('assigned', 0):,} assigned reads, "
                f"{counts_info.get('genes_detected', 0):,} genes detected."
            )

        data = {
            "output_dir": output_dir,
            "steps": step_results,
        }

        return ToolResult(
            success=True,
            data=data,
            summary=" ".join(summary_parts),
            files=output_files,
            warnings=warnings,
        )

    def _run_fastqc(
        self, reads_1: str, reads_2: str | None, qc_dir: str
    ) -> dict:
        """Run FastQC on input reads and return parsed results."""
        cmd = ["fastqc", reads_1, "-o", qc_dir, "--extract", "--quiet"]
        if reads_2:
            cmd.insert(2, reads_2)

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        output: dict = {"status": "success" if result.returncode == 0 else "failed"}

        if result.returncode != 0:
            output["error"] = result.stderr
            return output

        # Parse summary for R1
        stem = Path(reads_1).name.replace(".fastq.gz", "").replace(".fastq", "")
        summary_file = Path(qc_dir) / f"{stem}_fastqc" / "summary.txt"

        modules: dict[str, str] = {}
        qc_warnings: list[str] = []
        if summary_file.exists():
            for line in summary_file.read_text().splitlines():
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    modules[parts[1]] = parts[0]
            failed = [m for m, s in modules.items() if s == "FAIL"]
            if failed:
                qc_warnings.append(f"Failed modules: {', '.join(failed)}")

        output["modules"] = modules
        output["warnings"] = qc_warnings
        return output

    def _run_hisat2(
        self,
        reads_1: str,
        reads_2: str | None,
        reference: str,
        threads: int,
        output_bam: str,
    ) -> dict:
        """Run HISAT2 alignment and return statistics."""
        # Build index if needed
        index_base = reference.rsplit(".", 1)[0] if "." in reference else reference
        hisat2_index_files = [f"{index_base}.1.ht2", f"{index_base}.1.ht2l"]
        index_exists = any(Path(f).exists() for f in hisat2_index_files)

        if not index_exists:
            logger.info("Building HISAT2 index...")
            idx_result = subprocess.run(
                ["hisat2-build", reference, index_base],
                capture_output=True, text=True, timeout=7200,
            )
            if idx_result.returncode != 0:
                return {"status": "failed", "error": f"Index build failed: {idx_result.stderr}"}

        # Run HISAT2
        hisat2_cmd = ["hisat2", "-x", index_base, "-p", str(threads)]
        if reads_2:
            hisat2_cmd.extend(["-1", reads_1, "-2", reads_2])
        else:
            hisat2_cmd.extend(["-U", reads_1])

        sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", output_bam, "-"]

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
            return {"status": "failed", "error": f"Sorting failed: {sort_stderr.decode()}"}

        # Index BAM
        subprocess.run(
            ["samtools", "index", output_bam],
            capture_output=True, text=True, timeout=600,
        )

        # Collect stats with flagstat
        flagstat = subprocess.run(
            ["samtools", "flagstat", output_bam],
            capture_output=True, text=True, timeout=300,
        )

        total_reads = 0
        mapped_reads = 0
        if flagstat.returncode == 0:
            for line in flagstat.stdout.splitlines():
                if "in total" in line:
                    total_reads = int(line.split()[0])
                elif "mapped (" in line and "primary" not in line:
                    mapped_reads = int(line.split()[0])

        mapping_rate = (mapped_reads / total_reads * 100) if total_reads > 0 else 0.0

        return {
            "status": "success",
            "output_bam": output_bam,
            "total_reads": total_reads,
            "mapped_reads": mapped_reads,
            "mapping_rate": round(mapping_rate, 2),
            "hisat2_log": hisat2_stderr[-2000:] if hisat2_stderr else "",
        }

    def _run_featurecounts(
        self,
        bam_file: str,
        annotation_gtf: str,
        counts_file: str,
        threads: int,
        paired: bool,
    ) -> dict:
        """Run featureCounts and return assignment statistics."""
        cmd = [
            "featureCounts",
            "-a", annotation_gtf,
            "-o", counts_file,
            "-T", str(threads),
            "-g", "gene_id",
        ]
        if paired:
            cmd.append("-p")
            cmd.append("--countReadPairs")
        cmd.append(bam_file)

        logger.info(f"Running: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        output: dict = {"status": "success" if result.returncode == 0 else "failed"}

        if result.returncode != 0:
            output["error"] = result.stderr
            return output

        output["counts_file"] = counts_file

        # Parse the featureCounts summary file
        summary_file = counts_file + ".summary"
        if Path(summary_file).exists():
            assigned = 0
            unassigned = 0
            for line in Path(summary_file).read_text().splitlines():
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[0] != "Status":
                    count = int(parts[1])
                    if parts[0] == "Assigned":
                        assigned = count
                    else:
                        unassigned += count
            output["assigned"] = assigned
            output["unassigned"] = unassigned
            total = assigned + unassigned
            output["assignment_rate"] = round(
                (assigned / total * 100) if total > 0 else 0.0, 2
            )

        # Count detected genes (non-zero counts)
        counts_path = Path(counts_file)
        if counts_path.exists():
            genes_detected = 0
            total_genes = 0
            for line in counts_path.read_text().splitlines():
                if line.startswith("#") or line.startswith("Geneid"):
                    continue
                total_genes += 1
                parts = line.split("\t")
                if parts and int(parts[-1]) > 0:
                    genes_detected += 1
            output["total_genes"] = total_genes
            output["genes_detected"] = genes_detected

        return output
