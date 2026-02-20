"""ChEMBL and PubChem compound search tool."""

from __future__ import annotations

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.drugdiscovery.compound_search")


class CompoundSearchTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="compound_search",
            description=(
                "Search for compounds in ChEMBL or PubChem databases by name, synonym, "
                "or identifier. Returns compound names, database IDs, SMILES strings, "
                "and basic molecular properties."
            ),
            parameters=[
                ToolParameter(
                    name="query",
                    type=ParameterType.STRING,
                    description="Search query: compound name, synonym, or identifier (e.g., 'aspirin', 'ibuprofen')",
                ),
                ToolParameter(
                    name="database",
                    type=ParameterType.STRING,
                    description="Database to search",
                    enum=["chembl", "pubchem"],
                ),
                ToolParameter(
                    name="limit",
                    type=ParameterType.INTEGER,
                    description="Maximum number of results to return",
                    required=False,
                    default=10,
                ),
            ],
            category="drugdiscovery",
        )

    def check_available(self) -> bool:
        try:
            import chembl_webresource_client  # noqa: F401
            import pubchempy  # noqa: F401

            return True
        except ImportError:
            return False

    async def _search_chembl(self, query: str, limit: int) -> ToolResult:
        """Search ChEMBL for compounds."""
        try:
            from chembl_webresource_client.new_client import new_client
        except ImportError:
            return ToolResult(
                success=False,
                error="chembl_webresource_client is not installed. Install with: pip install chembl_webresource_client",
                summary="ChEMBL client not available",
            )

        molecule = new_client.molecule
        results = molecule.search(query)

        compounds = []
        for i, mol in enumerate(results):
            if i >= limit:
                break
            struct = mol.get("molecule_structures") or {}
            props = mol.get("molecule_properties") or {}
            compounds.append({
                "name": mol.get("pref_name") or "N/A",
                "chembl_id": mol.get("molecule_chembl_id", ""),
                "smiles": struct.get("canonical_smiles", ""),
                "formula": props.get("full_molformula", ""),
                "molecular_weight": props.get("full_mwt"),
                "logp": props.get("alogp"),
                "max_phase": mol.get("max_phase"),
                "molecule_type": mol.get("molecule_type", ""),
            })

        if not compounds:
            return ToolResult(
                success=True,
                data={"compounds": [], "database": "chembl", "query": query},
                summary=f"No compounds found in ChEMBL for query '{query}'.",
            )

        summary_lines = [f"Found {len(compounds)} compound(s) in ChEMBL for '{query}':"]
        for c in compounds[:5]:
            name = c["name"] if c["name"] != "N/A" else c["chembl_id"]
            summary_lines.append(f"  - {name} (MW={c['molecular_weight']}, Phase={c['max_phase']})")
        if len(compounds) > 5:
            summary_lines.append(f"  ... and {len(compounds) - 5} more")

        return ToolResult(
            success=True,
            data={"compounds": compounds, "database": "chembl", "query": query, "count": len(compounds)},
            summary="\n".join(summary_lines),
        )

    async def _search_pubchem(self, query: str, limit: int) -> ToolResult:
        """Search PubChem for compounds."""
        try:
            import pubchempy as pcp
        except ImportError:
            return ToolResult(
                success=False,
                error="pubchempy is not installed. Install with: pip install pubchempy",
                summary="PubChemPy not available",
            )

        results = pcp.get_compounds(query, "name", listkey_count=limit)

        compounds = []
        for comp in results[:limit]:
            compounds.append({
                "name": comp.iupac_name or "N/A",
                "cid": comp.cid,
                "smiles": comp.isomeric_smiles or comp.canonical_smiles or "",
                "formula": comp.molecular_formula or "",
                "molecular_weight": comp.molecular_weight,
                "logp": comp.xlogp,
                "exact_mass": comp.exact_mass,
                "synonyms": (comp.synonyms or [])[:5],
            })

        if not compounds:
            return ToolResult(
                success=True,
                data={"compounds": [], "database": "pubchem", "query": query},
                summary=f"No compounds found in PubChem for query '{query}'.",
            )

        summary_lines = [f"Found {len(compounds)} compound(s) in PubChem for '{query}':"]
        for c in compounds[:5]:
            name = c["name"] if c["name"] != "N/A" else f"CID {c['cid']}"
            summary_lines.append(f"  - {name} (CID={c['cid']}, MW={c['molecular_weight']})")
        if len(compounds) > 5:
            summary_lines.append(f"  ... and {len(compounds) - 5} more")

        return ToolResult(
            success=True,
            data={"compounds": compounds, "database": "pubchem", "query": query, "count": len(compounds)},
            summary="\n".join(summary_lines),
        )

    async def execute(self, **kwargs) -> ToolResult:
        query: str = kwargs["query"]
        database: str = kwargs["database"]
        limit: int = kwargs.get("limit", 10)

        if database == "chembl":
            return await self._search_chembl(query, limit)
        elif database == "pubchem":
            return await self._search_pubchem(query, limit)
        else:
            return ToolResult(
                success=False,
                error=f"Unknown database: {database}. Use 'chembl' or 'pubchem'.",
                summary=f"Invalid database selection: {database}",
            )
