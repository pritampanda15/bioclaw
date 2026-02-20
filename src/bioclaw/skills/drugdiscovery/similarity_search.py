"""Molecular fingerprint similarity search using RDKit."""

from __future__ import annotations

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.drugdiscovery.similarity_search")


class SimilaritySearchTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="similarity_search",
            description=(
                "Compute molecular fingerprint similarity between a query molecule and a "
                "set of reference molecules. Uses Tanimoto coefficient on ECFP4 (Morgan) or "
                "MACCS key fingerprints. Returns ranked similarity scores and identifies "
                "compounds above the similarity threshold."
            ),
            parameters=[
                ToolParameter(
                    name="smiles",
                    type=ParameterType.STRING,
                    description="SMILES string of the query molecule",
                ),
                ToolParameter(
                    name="reference_smiles",
                    type=ParameterType.ARRAY,
                    description="List of SMILES strings to compare against",
                    items_type=ParameterType.STRING,
                ),
                ToolParameter(
                    name="fingerprint_type",
                    type=ParameterType.STRING,
                    description="Type of molecular fingerprint to use",
                    required=False,
                    default="ecfp4",
                    enum=["ecfp4", "maccs"],
                ),
                ToolParameter(
                    name="threshold",
                    type=ParameterType.NUMBER,
                    description="Minimum Tanimoto similarity score to report (0.0-1.0)",
                    required=False,
                    default=0.7,
                ),
            ],
            category="drugdiscovery",
        )

    def check_available(self) -> bool:
        try:
            from rdkit import Chem  # noqa: F401
            from rdkit.Chem import AllChem, MACCSkeys  # noqa: F401

            return True
        except ImportError:
            return False

    async def execute(self, **kwargs) -> ToolResult:
        smiles: str = kwargs["smiles"]
        reference_smiles: list[str] = kwargs["reference_smiles"]
        fingerprint_type: str = kwargs.get("fingerprint_type", "ecfp4")
        threshold: float = kwargs.get("threshold", 0.7)

        try:
            from rdkit import Chem, DataStructs
            from rdkit.Chem import AllChem, MACCSkeys
        except ImportError:
            return ToolResult(
                success=False,
                error="RDKit is not installed. Install with: pip install rdkit",
                summary="RDKit not available",
            )

        # Parse query molecule
        query_mol = Chem.MolFromSmiles(smiles)
        if query_mol is None:
            return ToolResult(
                success=False,
                error=f"Invalid query SMILES: {smiles}",
                summary=f"Could not parse query SMILES: {smiles}",
            )

        # Generate query fingerprint
        if fingerprint_type == "ecfp4":
            query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, radius=2, nBits=2048)
        else:
            query_fp = MACCSkeys.GenMACCSKeys(query_mol)

        # Compare against references
        results = []
        warnings: list[str] = []
        for i, ref_smi in enumerate(reference_smiles):
            ref_mol = Chem.MolFromSmiles(ref_smi)
            if ref_mol is None:
                warnings.append(f"Could not parse reference SMILES at index {i}: {ref_smi}")
                continue

            if fingerprint_type == "ecfp4":
                ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, radius=2, nBits=2048)
            else:
                ref_fp = MACCSkeys.GenMACCSKeys(ref_mol)

            similarity = round(DataStructs.TanimotoSimilarity(query_fp, ref_fp), 4)
            results.append({
                "reference_smiles": ref_smi,
                "similarity": similarity,
                "above_threshold": similarity >= threshold,
                "index": i,
            })

        # Sort by similarity descending
        results.sort(key=lambda x: x["similarity"], reverse=True)

        similar = [r for r in results if r["above_threshold"]]

        data = {
            "query_smiles": smiles,
            "fingerprint_type": fingerprint_type,
            "threshold": threshold,
            "total_references": len(reference_smiles),
            "parsed_references": len(results),
            "above_threshold_count": len(similar),
            "results": results,
        }

        # Build summary
        summary_lines = [
            f"Similarity search ({fingerprint_type.upper()}, threshold={threshold}): "
            f"{len(similar)}/{len(results)} compounds above threshold."
        ]
        for r in results[:5]:
            marker = "*" if r["above_threshold"] else " "
            summary_lines.append(f"  {marker} Index {r['index']}: similarity={r['similarity']}")
        if len(results) > 5:
            summary_lines.append(f"  ... and {len(results) - 5} more")

        return ToolResult(
            success=True,
            data=data,
            summary="\n".join(summary_lines),
            warnings=warnings,
        )
