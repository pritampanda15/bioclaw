"""Molecular property calculation tool using RDKit."""

from __future__ import annotations

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult


class MolecularPropertiesTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="molecular_properties",
            description=(
                "Calculate molecular properties for a compound given its SMILES string. "
                "Returns molecular weight, LogP, hydrogen bond donors/acceptors, TPSA, "
                "rotatable bonds, Lipinski rule-of-five assessment, QED drug-likeness score, "
                "and PAINS filter results."
            ),
            parameters=[
                ToolParameter(
                    name="smiles",
                    type=ParameterType.STRING,
                    description="SMILES string of the molecule (e.g., 'CC(=O)Oc1ccccc1C(=O)O' for aspirin)",
                ),
                ToolParameter(
                    name="include_pains",
                    type=ParameterType.BOOLEAN,
                    description="Whether to run PAINS filter check",
                    required=False,
                    default=True,
                ),
            ],
            category="drugdiscovery",
        )

    def check_available(self) -> bool:
        try:
            from rdkit import Chem  # noqa: F401

            return True
        except ImportError:
            return False

    async def execute(self, **kwargs) -> ToolResult:
        smiles: str = kwargs["smiles"]
        include_pains: bool = kwargs.get("include_pains", True)

        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, FilterCatalog, QED, rdMolDescriptors
        except ImportError:
            return ToolResult(
                success=False,
                error="RDKit is not installed. Install with: pip install rdkit",
                summary="RDKit not available",
            )

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ToolResult(
                success=False,
                error=f"Invalid SMILES: {smiles}",
                summary=f"Could not parse SMILES string: {smiles}",
            )

        # Basic descriptors
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        tpsa = round(Descriptors.TPSA(mol), 2)
        rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        formula = rdMolDescriptors.CalcMolFormula(mol)

        # Lipinski Rule of Five
        lipinski_violations = sum([
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10,
        ])
        lipinski_pass = lipinski_violations <= 1

        # QED drug-likeness
        qed_score = round(QED.qed(mol), 3)

        # PAINS filter
        pains_alerts = []
        if include_pains:
            params = FilterCatalog.FilterCatalogParams()
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
            catalog = FilterCatalog.FilterCatalog(params)
            entry = catalog.GetFirstMatch(mol)
            if entry is not None:
                pains_alerts.append(entry.GetDescription())

        data = {
            "smiles": smiles,
            "formula": formula,
            "molecular_weight": mw,
            "logp": logp,
            "hbd": hbd,
            "hba": hba,
            "tpsa": tpsa,
            "rotatable_bonds": rotatable_bonds,
            "num_rings": num_rings,
            "num_aromatic_rings": num_aromatic_rings,
            "heavy_atoms": heavy_atoms,
            "lipinski": {
                "pass": lipinski_pass,
                "violations": lipinski_violations,
            },
            "qed": qed_score,
            "pains_alerts": pains_alerts,
        }

        warnings = []
        if not lipinski_pass:
            warnings.append(f"Lipinski Ro5: {lipinski_violations} violations")
        if pains_alerts:
            warnings.append(f"PAINS alerts: {', '.join(pains_alerts)}")

        summary = (
            f"Properties of {formula} (MW={mw}, LogP={logp}, HBD={hbd}, HBA={hba}, "
            f"TPSA={tpsa}, QED={qed_score}). "
            f"Lipinski: {'PASS' if lipinski_pass else f'FAIL ({lipinski_violations} violations)'}. "
            f"PAINS: {'clean' if not pains_alerts else ', '.join(pains_alerts)}."
        )

        return ToolResult(success=True, data=data, summary=summary, warnings=warnings)
