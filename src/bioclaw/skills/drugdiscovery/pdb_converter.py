"""PDB to PDBQT converter tool for receptor preparation."""

from __future__ import annotations

from pathlib import Path

from bioclaw.config.settings import get_settings
from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.drugdiscovery.pdb_converter")

# Common ions found in PDB crystal structures
COMMON_IONS = {
    "NA", "CL", "MG", "ZN", "CA", "FE", "MN", "CO", "NI", "CU",
    "K", "BR", "IOD", "CD", "HG", "SR", "BA", "CS", "LI", "RB",
    "SO4", "PO4", "NO3",
}

# AutoDock atom type mapping for receptor atoms.
# Vina computes its own scoring internally, so charges are set to 0.000.
# The atom types guide which interaction terms Vina applies.
_AD_ATOM_TYPES = {
    "C": "C",
    "N": "N",
    "O": "OA",
    "S": "SA",
    "H": "HD",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
    "SE": "Se",
    "FE": "Fe",
    "ZN": "Zn",
    "MG": "Mg",
    "MN": "Mn",
    "CA": "Ca",
}


def clean_pdb(
    pdb_path: Path,
    *,
    remove_waters: bool = True,
    remove_ligands: bool = True,
    remove_ions: bool = True,
) -> str:
    """Read a PDB file and return cleaned PDB content as a string.

    Filters out waters, ions, and non-protein HETATM records based on flags.
    Keeps ATOM records (protein backbone/side-chains) unconditionally.
    """
    kept_lines: list[str] = []

    for line in pdb_path.read_text().splitlines():
        record = line[:6].strip()

        if record == "ATOM":
            kept_lines.append(line)
        elif record == "HETATM":
            resname = line[17:20].strip()

            if remove_waters and resname == "HOH":
                continue
            if remove_ions and resname in COMMON_IONS:
                continue
            if remove_ligands:
                continue

            kept_lines.append(line)
        elif record in ("TER", "END", "MODEL", "ENDMDL", "HEADER", "TITLE",
                        "COMPND", "REMARK", "CRYST1", "SCALE1", "SCALE2",
                        "SCALE3", "ORIGX1", "ORIGX2", "ORIGX3"):
            kept_lines.append(line)

    kept_lines.append("END")
    return "\n".join(kept_lines) + "\n"


def _element_from_pdb_line(line: str) -> str:
    """Extract the element symbol from a PDB ATOM/HETATM line.

    Uses columns 77-78 (element field) when present, otherwise infers
    from the atom name in columns 13-16.
    """
    if len(line) >= 78:
        elem = line[76:78].strip()
        if elem:
            return elem.upper()

    # Fallback: derive from atom name (cols 13-16).
    # Standard PDB atom names: column 13 is blank for 1-3 char names,
    # or the first char of 4-char names (e.g. " CA " vs "1HG ").
    atom_name = line[12:16].strip()
    # Strip leading digits (e.g. "1HG" → "HG")
    stripped = atom_name.lstrip("0123456789")
    if not stripped:
        return "C"  # shouldn't happen, safe fallback

    # Two-char elements that appear in proteins: FE, ZN, MG, etc.
    if len(stripped) >= 2 and stripped[:2] in _AD_ATOM_TYPES:
        return stripped[:2]
    return stripped[0].upper()


def _get_ad_type(element: str) -> str:
    """Map an element symbol to its AutoDock atom type."""
    return _AD_ATOM_TYPES.get(element, element[:2])


def pdb_string_to_pdbqt(pdb_string: str) -> str:
    """Convert a cleaned PDB string to PDBQT format for receptor use.

    Assigns AutoDock atom types based on element and writes PDBQT lines
    with zero partial charges (Vina computes its own scoring).
    This is a receptor-oriented converter — Meeko's MoleculePreparation
    is designed for ligands and does not handle proteins correctly.
    """
    output_lines: list[str] = []

    for line in pdb_string.splitlines():
        record = line[:6].strip()

        if record in ("ATOM", "HETATM"):
            element = _element_from_pdb_line(line)
            ad_type = _get_ad_type(element)

            # Pad base line to 54 chars (coordinates end)
            base = line[:54].ljust(54)
            # Occupancy + B-factor (cols 55-66), keep original or default
            occ_bfac = line[54:66] if len(line) >= 66 else "  1.00  0.00"

            pdbqt_line = f"{base}{occ_bfac}    {0.000:+.3f} {ad_type:<2s}"
            output_lines.append(pdbqt_line)
        elif record in ("TER", "END", "MODEL", "ENDMDL", "REMARK"):
            output_lines.append(line)

    if not output_lines or output_lines[-1].strip() != "END":
        output_lines.append("END")

    return "\n".join(output_lines) + "\n"


class PDBConverterTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="pdb_to_pdbqt",
            description=(
                "Convert a PDB file to PDBQT format for molecular docking. "
                "Cleans the structure by removing waters, ions, and non-protein "
                "ligands (configurable), then assigns AutoDock atom types and "
                "writes PDBQT format suitable for receptor use with Vina."
            ),
            parameters=[
                ToolParameter(
                    name="pdb_file",
                    type=ParameterType.STRING,
                    description="Path to the input PDB file",
                ),
                ToolParameter(
                    name="remove_waters",
                    type=ParameterType.BOOLEAN,
                    description="Remove water molecules (HOH)",
                    required=False,
                    default=True,
                ),
                ToolParameter(
                    name="remove_ligands",
                    type=ParameterType.BOOLEAN,
                    description="Remove non-protein ligands (HETATM records)",
                    required=False,
                    default=True,
                ),
                ToolParameter(
                    name="remove_ions",
                    type=ParameterType.BOOLEAN,
                    description="Remove common ions (NA, CL, MG, ZN, etc.)",
                    required=False,
                    default=True,
                ),
                ToolParameter(
                    name="output_file",
                    type=ParameterType.STRING,
                    description="Path for the output PDBQT file (default: bioclaw_data/structures/<name>_clean.pdbqt)",
                    required=False,
                ),
            ],
            category="drugdiscovery",
        )

    def check_available(self) -> bool:
        # Pure Python PDB text processing — no external dependencies needed.
        return True

    async def execute(self, **kwargs) -> ToolResult:
        pdb_file: str = kwargs["pdb_file"]
        remove_waters: bool = kwargs.get("remove_waters", True)
        remove_ligands: bool = kwargs.get("remove_ligands", True)
        remove_ions: bool = kwargs.get("remove_ions", True)
        output_file: str | None = kwargs.get("output_file")

        pdb_path = Path(pdb_file)
        if not pdb_path.exists():
            return ToolResult(
                success=False,
                error=f"PDB file not found: {pdb_file}",
                summary=f"PDB file does not exist: {pdb_file}",
            )

        # Determine output path
        if output_file:
            output_path = Path(output_file)
        else:
            settings = get_settings()
            structures_dir = settings.data_dir / "structures"
            structures_dir.mkdir(parents=True, exist_ok=True)
            output_path = structures_dir / f"{pdb_path.stem}_clean.pdbqt"

        try:
            # Clean the PDB
            cleaned_pdb = clean_pdb(
                pdb_path,
                remove_waters=remove_waters,
                remove_ligands=remove_ligands,
                remove_ions=remove_ions,
            )

            # Convert to PDBQT
            pdbqt_string = pdb_string_to_pdbqt(cleaned_pdb)

            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(pdbqt_string)

            return ToolResult(
                success=True,
                data={
                    "input_file": str(pdb_path),
                    "output_file": str(output_path),
                    "remove_waters": remove_waters,
                    "remove_ligands": remove_ligands,
                    "remove_ions": remove_ions,
                },
                summary=f"Converted {pdb_path.name} to PDBQT: {output_path}",
                files=[output_path],
            )

        except Exception as e:
            return ToolResult(
                success=False,
                error=f"PDB to PDBQT conversion failed: {e}",
                summary=f"Conversion failed: {e}",
            )
