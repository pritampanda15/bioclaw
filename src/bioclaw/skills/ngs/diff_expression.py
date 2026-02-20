"""Differential expression analysis tool using PyDESeq2."""

from __future__ import annotations

from pathlib import Path

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.ngs.diff_expression")


class DiffExpressionTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="diff_expression",
            description=(
                "Perform differential expression analysis on an RNA-seq count matrix "
                "using PyDESeq2. Takes a count matrix (genes x samples) and sample "
                "metadata file, runs DESeq2 normalization and statistical testing, "
                "and returns significant differentially expressed genes with fold "
                "changes and adjusted p-values."
            ),
            parameters=[
                ToolParameter(
                    name="counts_file",
                    type=ParameterType.STRING,
                    description=(
                        "Path to the count matrix file (CSV or TSV). Rows are genes, "
                        "columns are samples. First column should be gene IDs."
                    ),
                ),
                ToolParameter(
                    name="metadata_file",
                    type=ParameterType.STRING,
                    description=(
                        "Path to the sample metadata file (CSV or TSV). Must contain "
                        "a sample ID column matching the count matrix column names."
                    ),
                ),
                ToolParameter(
                    name="condition_column",
                    type=ParameterType.STRING,
                    description="Column name in metadata that defines the experimental conditions",
                ),
                ToolParameter(
                    name="reference_condition",
                    type=ParameterType.STRING,
                    description="The reference/control condition for comparison",
                ),
                ToolParameter(
                    name="alpha",
                    type=ParameterType.NUMBER,
                    description="Significance threshold for adjusted p-value",
                    required=False,
                    default=0.05,
                ),
            ],
            category="ngs",
            is_long_running=True,
            required_binaries=[],
        )

    def check_available(self) -> bool:
        try:
            import pandas  # noqa: F401
            import pydeseq2  # noqa: F401

            return True
        except ImportError:
            return False

    async def execute(self, **kwargs) -> ToolResult:
        counts_file = kwargs["counts_file"]
        metadata_file = kwargs["metadata_file"]
        condition_column = kwargs["condition_column"]
        reference_condition = kwargs["reference_condition"]
        alpha = kwargs.get("alpha", 0.05)

        # Validate input files
        for label, path_str in [("Counts file", counts_file), ("Metadata file", metadata_file)]:
            if not Path(path_str).exists():
                return ToolResult(
                    success=False,
                    error=f"{label} not found: {path_str}",
                    summary=f"{label} does not exist: {path_str}",
                )

        # Import required libraries
        try:
            import pandas as pd
        except ImportError:
            return ToolResult(
                success=False,
                error="pandas is not installed. Install with: pip install pandas",
                summary="pandas not available",
            )

        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError:
            return ToolResult(
                success=False,
                error="PyDESeq2 is not installed. Install with: pip install pydeseq2",
                summary="PyDESeq2 not available",
            )

        # Load data
        try:
            counts_path = Path(counts_file)
            sep = "\t" if counts_path.suffix in (".tsv", ".txt") else ","
            counts_df = pd.read_csv(counts_file, sep=sep, index_col=0)

            meta_path = Path(metadata_file)
            meta_sep = "\t" if meta_path.suffix in (".tsv", ".txt") else ","
            metadata_df = pd.read_csv(metadata_file, sep=meta_sep, index_col=0)
        except Exception as e:
            return ToolResult(
                success=False,
                error=f"Failed to load input files: {e}",
                summary=f"Data loading error: {e}",
            )

        # Validate condition column
        if condition_column not in metadata_df.columns:
            available = ", ".join(metadata_df.columns.tolist())
            return ToolResult(
                success=False,
                error=f"Column '{condition_column}' not found in metadata. Available: {available}",
                summary=f"Condition column '{condition_column}' not in metadata",
            )

        # Validate reference condition
        conditions = metadata_df[condition_column].unique().tolist()
        if reference_condition not in conditions:
            return ToolResult(
                success=False,
                error=(
                    f"Reference condition '{reference_condition}' not found. "
                    f"Available conditions: {conditions}"
                ),
                summary=f"Reference condition '{reference_condition}' not in metadata",
            )

        # Ensure sample overlap
        common_samples = counts_df.columns.intersection(metadata_df.index)
        if len(common_samples) == 0:
            return ToolResult(
                success=False,
                error="No common samples between count matrix columns and metadata index",
                summary="No overlapping sample IDs between counts and metadata",
            )

        counts_df = counts_df[common_samples]
        metadata_df = metadata_df.loc[common_samples]

        # Transpose counts: PyDESeq2 expects samples x genes
        counts_transposed = counts_df.T
        # Ensure integer counts
        counts_transposed = counts_transposed.round().astype(int)

        warnings = []
        if len(common_samples) < len(counts_df.columns):
            dropped = len(counts_df.columns) - len(common_samples)
            warnings.append(f"{dropped} samples in count matrix not found in metadata")

        try:
            # Run DESeq2
            dds = DeseqDataSet(
                counts=counts_transposed,
                metadata=metadata_df,
                design_factors=condition_column,
                refit_cooks=True,
                ref_level=[condition_column, reference_condition],
            )
            dds.deseq2()

            # Statistical results
            stat_res = DeseqStats(dds, alpha=alpha)
            stat_res.summary()
            results_df = stat_res.results_df

        except Exception as e:
            return ToolResult(
                success=False,
                error=f"DESeq2 analysis failed: {e}",
                summary=f"DESeq2 execution error: {e}",
            )

        # Process results
        try:
            results_df = results_df.dropna(subset=["padj"])

            sig_genes = results_df[results_df["padj"] < alpha].copy()
            sig_up = sig_genes[sig_genes["log2FoldChange"] > 0]
            sig_down = sig_genes[sig_genes["log2FoldChange"] < 0]

            # Top significant genes by adjusted p-value
            top_genes = sig_genes.nsmallest(20, "padj")
            top_genes_list = []
            for gene_id, row in top_genes.iterrows():
                top_genes_list.append({
                    "gene_id": str(gene_id),
                    "log2_fold_change": round(float(row["log2FoldChange"]), 3),
                    "p_adjusted": float(row["padj"]),
                    "base_mean": round(float(row["baseMean"]), 2),
                })

            data = {
                "counts_file": counts_file,
                "metadata_file": metadata_file,
                "condition_column": condition_column,
                "reference_condition": reference_condition,
                "alpha": alpha,
                "total_genes_tested": len(results_df),
                "significant_genes": len(sig_genes),
                "upregulated": len(sig_up),
                "downregulated": len(sig_down),
                "conditions_compared": [
                    c for c in conditions if c != reference_condition
                ],
                "top_significant_genes": top_genes_list,
            }

            summary = (
                f"DE analysis complete: {len(results_df):,} genes tested, "
                f"{len(sig_genes):,} significant at alpha={alpha} "
                f"({len(sig_up):,} up, {len(sig_down):,} down). "
                f"Reference: {reference_condition}."
            )

            if len(sig_genes) == 0:
                warnings.append(
                    "No significant DE genes found. Consider adjusting alpha or "
                    "checking sample quality/replicates."
                )

            return ToolResult(
                success=True,
                data=data,
                summary=summary,
                warnings=warnings,
            )

        except Exception as e:
            return ToolResult(
                success=False,
                error=f"Results processing failed: {e}",
                summary=f"Failed to process DE results: {e}",
            )
