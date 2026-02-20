"""Murcko scaffold analysis and R-group decomposition tool."""

from __future__ import annotations

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.drugdiscovery.scaffold_analysis")


class ScaffoldAnalysisTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="scaffold_analysis",
            description=(
                "Analyze molecular scaffolds for a set of compounds. Murcko mode extracts "
                "Bemis-Murcko scaffolds and generic frameworks, identifying shared scaffolds "
                "across the compound set. R-group mode performs R-group decomposition around "
                "a provided core structure, identifying substituent variations."
            ),
            parameters=[
                ToolParameter(
                    name="smiles_list",
                    type=ParameterType.ARRAY,
                    description="List of SMILES strings to analyze",
                    items_type=ParameterType.STRING,
                ),
                ToolParameter(
                    name="analysis_type",
                    type=ParameterType.STRING,
                    description="Type of scaffold analysis to perform",
                    required=False,
                    default="murcko",
                    enum=["murcko", "rgroup"],
                ),
                ToolParameter(
                    name="core_smiles",
                    type=ParameterType.STRING,
                    description="Core structure SMILES for R-group decomposition (required for rgroup mode)",
                    required=False,
                ),
            ],
            category="drugdiscovery",
        )

    def check_available(self) -> bool:
        try:
            from rdkit import Chem  # noqa: F401
            from rdkit.Chem.Scaffolds import MurckoScaffold  # noqa: F401

            return True
        except ImportError:
            return False

    def _murcko_analysis(self, smiles_list: list[str]) -> ToolResult:
        """Extract Murcko scaffolds and generic frameworks."""
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        molecules = []
        warnings: list[str] = []
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                warnings.append(f"Could not parse SMILES at index {i}: {smi}")
            else:
                molecules.append((i, smi, mol))

        if not molecules:
            return ToolResult(
                success=False,
                error="No valid molecules provided",
                summary="Could not parse any of the provided SMILES strings.",
                warnings=warnings,
            )

        # Extract scaffolds
        scaffold_map: dict[str, list[dict]] = {}
        compound_results = []

        for idx, smi, mol in molecules:
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smi = Chem.MolToSmiles(scaffold)

                generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
                generic_smi = Chem.MolToSmiles(generic)
            except Exception as e:
                warnings.append(f"Scaffold extraction failed for index {idx}: {e}")
                scaffold_smi = ""
                generic_smi = ""

            entry = {
                "index": idx,
                "smiles": smi,
                "scaffold": scaffold_smi,
                "generic_framework": generic_smi,
            }
            compound_results.append(entry)

            if scaffold_smi:
                scaffold_map.setdefault(scaffold_smi, []).append({"index": idx, "smiles": smi})

        # Identify shared scaffolds (appearing in more than one molecule)
        shared_scaffolds = [
            {"scaffold": smi, "count": len(members), "members": members}
            for smi, members in scaffold_map.items()
            if len(members) > 1
        ]
        shared_scaffolds.sort(key=lambda x: x["count"], reverse=True)

        unique_scaffolds = len(scaffold_map)

        data = {
            "analysis_type": "murcko",
            "total_molecules": len(smiles_list),
            "parsed_molecules": len(molecules),
            "unique_scaffolds": unique_scaffolds,
            "shared_scaffolds": shared_scaffolds,
            "compounds": compound_results,
        }

        summary_lines = [
            f"Murcko scaffold analysis: {len(molecules)} molecules parsed, "
            f"{unique_scaffolds} unique scaffolds found.",
        ]
        if shared_scaffolds:
            summary_lines.append(f"Shared scaffolds ({len(shared_scaffolds)}):")
            for s in shared_scaffolds[:5]:
                summary_lines.append(f"  - {s['scaffold']} (n={s['count']})")
        else:
            summary_lines.append("No shared scaffolds found.")

        return ToolResult(
            success=True,
            data=data,
            summary="\n".join(summary_lines),
            warnings=warnings,
        )

    def _rgroup_analysis(self, smiles_list: list[str], core_smiles: str) -> ToolResult:
        """Perform R-group decomposition around a core structure."""
        from rdkit import Chem
        from rdkit.Chem import rdRGroupDecomposition

        core = Chem.MolFromSmiles(core_smiles)
        if core is None:
            return ToolResult(
                success=False,
                error=f"Invalid core SMILES: {core_smiles}",
                summary=f"Could not parse core SMILES: {core_smiles}",
            )

        molecules = []
        valid_smiles = []
        warnings: list[str] = []
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                warnings.append(f"Could not parse SMILES at index {i}: {smi}")
            else:
                molecules.append(mol)
                valid_smiles.append(smi)

        if not molecules:
            return ToolResult(
                success=False,
                error="No valid molecules provided",
                summary="Could not parse any of the provided SMILES strings.",
                warnings=warnings,
            )

        # Perform R-group decomposition
        try:
            rg_result, unmatched = rdRGroupDecomposition.RGroupDecompose(
                [core], molecules, asSmiles=True
            )
        except Exception as e:
            return ToolResult(
                success=False,
                error=f"R-group decomposition failed: {e}",
                summary=f"R-group decomposition failed: {e}",
            )

        if not rg_result:
            return ToolResult(
                success=True,
                data={
                    "analysis_type": "rgroup",
                    "core": core_smiles,
                    "matched": 0,
                    "unmatched": len(molecules),
                    "decompositions": [],
                },
                summary=f"No molecules matched the core structure: {core_smiles}",
                warnings=warnings,
            )

        # Build decomposition results
        decompositions = []
        r_group_keys = set()
        for i, entry in enumerate(rg_result):
            decomp = {"smiles": valid_smiles[i]}
            decomp.update(entry)
            decompositions.append(decomp)
            r_group_keys.update(k for k in entry.keys() if k.startswith("R"))

        # Analyze R-group diversity
        r_group_diversity = {}
        for key in sorted(r_group_keys):
            unique_groups = set(d.get(key, "") for d in decompositions if d.get(key))
            r_group_diversity[key] = {
                "unique_count": len(unique_groups),
                "groups": sorted(unique_groups),
            }

        matched_count = len(rg_result)
        unmatched_count = len(molecules) - matched_count

        data = {
            "analysis_type": "rgroup",
            "core": core_smiles,
            "matched": matched_count,
            "unmatched": unmatched_count,
            "r_groups": sorted(r_group_keys),
            "r_group_diversity": r_group_diversity,
            "decompositions": decompositions,
        }

        summary_lines = [
            f"R-group decomposition: {matched_count}/{len(molecules)} molecules matched "
            f"core '{core_smiles}'.",
            f"R-group positions: {', '.join(sorted(r_group_keys))}.",
        ]
        for key in sorted(r_group_keys):
            div = r_group_diversity[key]
            summary_lines.append(f"  {key}: {div['unique_count']} unique substituents")

        if unmatched_count > 0:
            warnings.append(f"{unmatched_count} molecule(s) did not match the core structure.")

        return ToolResult(
            success=True,
            data=data,
            summary="\n".join(summary_lines),
            warnings=warnings,
        )

    async def execute(self, **kwargs) -> ToolResult:
        smiles_list: list[str] = kwargs["smiles_list"]
        analysis_type: str = kwargs.get("analysis_type", "murcko")
        core_smiles: str | None = kwargs.get("core_smiles")

        if analysis_type == "rgroup":
            if not core_smiles:
                return ToolResult(
                    success=False,
                    error="core_smiles is required for R-group decomposition",
                    summary="Missing required parameter: core_smiles for rgroup analysis.",
                )
            return self._rgroup_analysis(smiles_list, core_smiles)
        else:
            return self._murcko_analysis(smiles_list)
