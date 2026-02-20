"""AutoDock Vina molecular docking wrapper tool."""

from __future__ import annotations

import asyncio
import json
import tempfile
from pathlib import Path

from bioclaw.config.settings import get_settings
from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.drugdiscovery.molecular_docking")


class MolecularDockingTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="molecular_docking",
            description=(
                "Perform molecular docking of a ligand into a receptor binding site using "
                "AutoDock Vina. Accepts a ligand SMILES string (converted to PDBQT via Meeko) "
                "and a receptor PDBQT file path. Returns binding affinity scores and docked poses."
            ),
            parameters=[
                ToolParameter(
                    name="ligand_smiles",
                    type=ParameterType.STRING,
                    description="SMILES string of the ligand molecule",
                ),
                ToolParameter(
                    name="receptor_pdbqt",
                    type=ParameterType.STRING,
                    description="Path to the receptor PDBQT file",
                ),
                ToolParameter(
                    name="center_x",
                    type=ParameterType.NUMBER,
                    description="X coordinate of the search box center (Angstroms)",
                ),
                ToolParameter(
                    name="center_y",
                    type=ParameterType.NUMBER,
                    description="Y coordinate of the search box center (Angstroms)",
                ),
                ToolParameter(
                    name="center_z",
                    type=ParameterType.NUMBER,
                    description="Z coordinate of the search box center (Angstroms)",
                ),
                ToolParameter(
                    name="size_x",
                    type=ParameterType.NUMBER,
                    description="Search box size in X dimension (Angstroms)",
                    required=False,
                    default=20,
                ),
                ToolParameter(
                    name="size_y",
                    type=ParameterType.NUMBER,
                    description="Search box size in Y dimension (Angstroms)",
                    required=False,
                    default=20,
                ),
                ToolParameter(
                    name="size_z",
                    type=ParameterType.NUMBER,
                    description="Search box size in Z dimension (Angstroms)",
                    required=False,
                    default=20,
                ),
                ToolParameter(
                    name="exhaustiveness",
                    type=ParameterType.INTEGER,
                    description="Exhaustiveness of the search (higher = more thorough but slower)",
                    required=False,
                    default=8,
                ),
            ],
            category="drugdiscovery",
            is_long_running=True,
            required_binaries=["vina"],
        )

    def check_available(self) -> bool:
        """Check for vina binary and required Python packages."""
        import shutil

        if shutil.which("vina") is None:
            return False
        try:
            import meeko  # noqa: F401
            from rdkit import Chem  # noqa: F401

            return True
        except ImportError:
            return False

    def _smiles_to_pdbqt(self, smiles: str, output_path: Path) -> Path:
        """Convert SMILES to PDBQT using RDKit and Meeko."""
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result != 0:
            AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)

        # Prepare with Meeko
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        # Write PDBQT
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            raise RuntimeError(f"Meeko PDBQT generation failed: {error_msg}")

        output_path.write_text(pdbqt_string)
        return output_path

    def _parse_vina_output(self, stdout: str) -> list[dict]:
        """Parse Vina output to extract binding modes."""
        modes = []
        in_table = False
        for line in stdout.splitlines():
            line = line.strip()
            if line.startswith("-----+"):
                in_table = True
                continue
            if in_table and line and line[0].isdigit():
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        modes.append({
                            "mode": int(parts[0]),
                            "affinity_kcal_mol": float(parts[1]),
                            "rmsd_lb": float(parts[2]),
                            "rmsd_ub": float(parts[3]),
                        })
                    except (ValueError, IndexError):
                        continue
            elif in_table and not line:
                break
        return modes

    async def execute(self, **kwargs) -> ToolResult:
        ligand_smiles: str = kwargs["ligand_smiles"]
        receptor_pdbqt: str = kwargs["receptor_pdbqt"]
        center_x: float = kwargs["center_x"]
        center_y: float = kwargs["center_y"]
        center_z: float = kwargs["center_z"]
        size_x: float = kwargs.get("size_x", 20)
        size_y: float = kwargs.get("size_y", 20)
        size_z: float = kwargs.get("size_z", 20)
        exhaustiveness: int = kwargs.get("exhaustiveness", 8)

        # Validate receptor file exists
        receptor_path = Path(receptor_pdbqt)
        if not receptor_path.exists():
            return ToolResult(
                success=False,
                error=f"Receptor file not found: {receptor_pdbqt}",
                summary=f"Receptor PDBQT file does not exist: {receptor_pdbqt}",
            )

        settings = get_settings()
        output_dir = settings.data_dir / "docking"
        output_dir.mkdir(parents=True, exist_ok=True)

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir_path = Path(tmpdir)

                # Convert SMILES to PDBQT
                ligand_pdbqt = tmpdir_path / "ligand.pdbqt"
                try:
                    self._smiles_to_pdbqt(ligand_smiles, ligand_pdbqt)
                except Exception as e:
                    return ToolResult(
                        success=False,
                        error=f"Ligand preparation failed: {e}",
                        summary=f"Could not convert SMILES to PDBQT: {e}",
                    )

                # Output file for docked poses
                output_pdbqt = output_dir / "docked_poses.pdbqt"

                # Build Vina command
                cmd = [
                    "vina",
                    "--receptor", str(receptor_path),
                    "--ligand", str(ligand_pdbqt),
                    "--center_x", str(center_x),
                    "--center_y", str(center_y),
                    "--center_z", str(center_z),
                    "--size_x", str(size_x),
                    "--size_y", str(size_y),
                    "--size_z", str(size_z),
                    "--exhaustiveness", str(exhaustiveness),
                    "--out", str(output_pdbqt),
                ]

                logger.info(f"Running Vina: {' '.join(cmd)}")
                proc = await asyncio.create_subprocess_exec(
                    *cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )
                stdout, stderr = await proc.communicate()

                stdout_text = stdout.decode("utf-8", errors="replace")
                stderr_text = stderr.decode("utf-8", errors="replace")

                if proc.returncode != 0:
                    return ToolResult(
                        success=False,
                        error=f"Vina exited with code {proc.returncode}: {stderr_text}",
                        summary=f"Docking failed: {stderr_text[:200]}",
                    )

            # Parse results
            modes = self._parse_vina_output(stdout_text)

            files = [output_pdbqt] if output_pdbqt.exists() else []

            data = {
                "ligand_smiles": ligand_smiles,
                "receptor": str(receptor_path),
                "search_box": {
                    "center": [center_x, center_y, center_z],
                    "size": [size_x, size_y, size_z],
                },
                "exhaustiveness": exhaustiveness,
                "modes": modes,
                "output_pdbqt": str(output_pdbqt) if output_pdbqt.exists() else None,
            }

            if modes:
                best = modes[0]
                summary = (
                    f"Docking completed: {len(modes)} binding modes found. "
                    f"Best affinity: {best['affinity_kcal_mol']} kcal/mol. "
                    f"Output saved to {output_pdbqt}."
                )
            else:
                summary = "Docking completed but no binding modes could be parsed from output."

            return ToolResult(
                success=True,
                data=data,
                summary=summary,
                files=files,
            )

        except Exception as e:
            return ToolResult(
                success=False,
                error=f"Docking failed: {e}",
                summary=f"Molecular docking encountered an error: {e}",
            )
