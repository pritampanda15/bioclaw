"""ADMET property prediction tool using RDKit descriptors."""

from __future__ import annotations

import math

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.drugdiscovery.admet_predictor")


class ADMETpredictorTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="admet_predictor",
            description=(
                "Predict ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) "
                "properties for a compound from its SMILES string. Uses RDKit molecular "
                "descriptors and rule-based models to estimate LogP, solubility (ESOL), "
                "TPSA, GI absorption, BBB permeability, CYP inhibition risk, and basic "
                "toxicity alerts."
            ),
            parameters=[
                ToolParameter(
                    name="smiles",
                    type=ParameterType.STRING,
                    description="SMILES string of the molecule to evaluate",
                ),
            ],
            category="drugdiscovery",
        )

    def check_available(self) -> bool:
        try:
            from rdkit import Chem  # noqa: F401
            from rdkit.Chem import Descriptors  # noqa: F401

            return True
        except ImportError:
            return False

    def _estimate_solubility(self, logp: float, mw: float, rotatable_bonds: int, ap: int) -> dict:
        """Estimate aqueous solubility using the ESOL model (Delaney 2004).

        log(S) = 0.16 - 0.63*logP - 0.0062*MW + 0.066*RB - 0.74*AP
        where AP = aromatic proportion (aromatic atoms / heavy atoms).
        """
        logs = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotatable_bonds - 0.74 * ap
        solubility_mg_ml = (10 ** logs) * mw

        if logs > -1:
            category = "highly_soluble"
        elif logs > -3:
            category = "soluble"
        elif logs > -5:
            category = "moderately_soluble"
        else:
            category = "poorly_soluble"

        return {
            "log_s": round(logs, 2),
            "solubility_mg_ml": round(solubility_mg_ml, 4),
            "category": category,
        }

    def _predict_bbb(self, tpsa: float, mw: float, logp: float, hbd: int) -> dict:
        """Predict BBB permeability using simple descriptor rules.

        Criteria: TPSA < 90, MW < 500, LogP 1-3, HBD <= 3.
        """
        score = 0
        reasons = []

        if tpsa < 90:
            score += 1
        else:
            reasons.append(f"TPSA={tpsa} > 90")

        if mw < 500:
            score += 1
        else:
            reasons.append(f"MW={mw} > 500")

        if 1.0 <= logp <= 3.0:
            score += 1
        else:
            reasons.append(f"LogP={logp} outside 1-3 range")

        if hbd <= 3:
            score += 1
        else:
            reasons.append(f"HBD={hbd} > 3")

        if score >= 3:
            prediction = "likely_permeable"
        elif score >= 2:
            prediction = "uncertain"
        else:
            prediction = "likely_impermeable"

        return {
            "prediction": prediction,
            "score": score,
            "max_score": 4,
            "limiting_factors": reasons,
        }

    def _predict_gi_absorption(self, tpsa: float, logp: float) -> dict:
        """Predict GI absorption using the BOILED-Egg model approximation.

        High absorption: TPSA <= 140 and -1 <= LogP <= 5.5.
        """
        high_absorption = tpsa <= 140 and -1 <= logp <= 5.5
        return {
            "prediction": "high" if high_absorption else "low",
            "tpsa": tpsa,
            "logp": logp,
        }

    def _predict_cyp_inhibition(self, logp: float, mw: float, hba: int, aromatic_rings: int) -> dict:
        """Predict CYP inhibition risk using simple descriptor rules.

        Higher risk: LogP > 3, MW > 300, multiple aromatic rings.
        """
        risk_score = 0
        flags = []

        if logp > 3.0:
            risk_score += 1
            flags.append("high_logp")
        if mw > 300:
            risk_score += 1
            flags.append("high_mw")
        if aromatic_rings >= 2:
            risk_score += 1
            flags.append("multiple_aromatic_rings")
        if hba > 5:
            risk_score += 1
            flags.append("high_hba")

        if risk_score >= 3:
            risk = "high"
        elif risk_score >= 2:
            risk = "moderate"
        else:
            risk = "low"

        return {
            "risk": risk,
            "risk_score": risk_score,
            "max_score": 4,
            "flags": flags,
        }

    def _assess_toxicity_alerts(self, mol) -> dict:
        """Check for common toxicity structural alerts using SMARTS patterns."""
        from rdkit import Chem

        alerts = {
            "nitro_aromatic": "[$(c[N+](=O)[O-]),$(c[N](=O)=O)]",
            "epoxide": "C1OC1",
            "michael_acceptor": "[C;!R]=[C;!R]-[C,N,O]=[O,N]",
            "acyl_halide": "[CX3](=[OX1])[F,Cl,Br,I]",
            "aldehyde": "[CX3H1](=O)[#6]",
            "alkyl_halide": "[CX4][F,Cl,Br,I]",
            "aniline": "c-[NX3;H2]",
            "hydrazine": "[NX3;H2]-[NX3;H2]",
        }

        found_alerts = []
        for name, smarts in alerts.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is not None and mol.HasSubstructMatch(pattern):
                found_alerts.append(name)

        return {
            "alerts_found": found_alerts,
            "alert_count": len(found_alerts),
            "classification": "flagged" if found_alerts else "clean",
        }

    async def execute(self, **kwargs) -> ToolResult:
        smiles: str = kwargs["smiles"]

        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors
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

        # Compute descriptors
        mw = round(Descriptors.MolWt(mol), 2)
        logp = round(Descriptors.MolLogP(mol), 2)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        tpsa = round(Descriptors.TPSA(mol), 2)
        rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_proportion = round(aromatic_atoms / heavy_atoms, 3) if heavy_atoms > 0 else 0

        # ADMET predictions
        absorption = {
            "gi_absorption": self._predict_gi_absorption(tpsa, logp),
            "solubility": self._estimate_solubility(logp, mw, rotatable_bonds, aromatic_proportion),
        }

        distribution = {
            "bbb_permeability": self._predict_bbb(tpsa, mw, logp, hbd),
            "tpsa": tpsa,
            "logp": logp,
        }

        metabolism = {
            "cyp_inhibition_risk": self._predict_cyp_inhibition(logp, mw, hba, aromatic_rings),
        }

        excretion = {
            "molecular_weight": mw,
            "mw_category": "low_clearance_risk" if mw < 500 else "high_clearance_risk",
            "rotatable_bonds": rotatable_bonds,
        }

        toxicity = self._assess_toxicity_alerts(mol)

        data = {
            "smiles": smiles,
            "descriptors": {
                "molecular_weight": mw,
                "logp": logp,
                "hbd": hbd,
                "hba": hba,
                "tpsa": tpsa,
                "rotatable_bonds": rotatable_bonds,
                "aromatic_rings": aromatic_rings,
                "heavy_atoms": heavy_atoms,
            },
            "absorption": absorption,
            "distribution": distribution,
            "metabolism": metabolism,
            "excretion": excretion,
            "toxicity": toxicity,
        }

        # Build summary
        gi = absorption["gi_absorption"]["prediction"]
        sol = absorption["solubility"]["category"]
        bbb = distribution["bbb_permeability"]["prediction"]
        cyp = metabolism["cyp_inhibition_risk"]["risk"]
        tox_count = toxicity["alert_count"]

        warnings = []
        if gi == "low":
            warnings.append("Low predicted GI absorption")
        if bbb == "likely_permeable":
            warnings.append("Likely BBB permeable (CNS exposure)")
        if cyp == "high":
            warnings.append("High CYP inhibition risk (drug-drug interactions)")
        if tox_count > 0:
            warnings.append(f"Toxicity alerts: {', '.join(toxicity['alerts_found'])}")
        if sol == "poorly_soluble":
            warnings.append("Poorly soluble compound")

        summary = (
            f"ADMET profile for {smiles}: "
            f"GI absorption={gi}, solubility={sol}, "
            f"BBB={bbb}, CYP inhibition risk={cyp}, "
            f"toxicity alerts={tox_count}. "
            f"Key descriptors: MW={mw}, LogP={logp}, TPSA={tpsa}, HBD={hbd}, HBA={hba}."
        )

        return ToolResult(
            success=True,
            data=data,
            summary=summary,
            warnings=warnings,
        )
