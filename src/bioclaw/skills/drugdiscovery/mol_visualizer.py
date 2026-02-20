"""2D/3D molecule rendering tool."""

from __future__ import annotations

import uuid
from pathlib import Path

from bioclaw.config.settings import get_settings
from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult
from bioclaw.utils.logging import get_logger

logger = get_logger("skills.drugdiscovery.mol_visualizer")


class MolVisualizerTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="mol_visualizer",
            description=(
                "Render a 2D or 3D visualization of a molecule from its SMILES string. "
                "2D mode produces a publication-quality structural diagram (PNG or SVG). "
                "3D mode generates an interactive HTML viewer using py3Dmol with a "
                "UFF-optimized 3D conformer."
            ),
            parameters=[
                ToolParameter(
                    name="smiles",
                    type=ParameterType.STRING,
                    description="SMILES string of the molecule to visualize",
                ),
                ToolParameter(
                    name="mode",
                    type=ParameterType.STRING,
                    description="Rendering mode: 2d for structural diagram, 3d for interactive viewer",
                    required=False,
                    default="2d",
                    enum=["2d", "3d"],
                ),
                ToolParameter(
                    name="output_format",
                    type=ParameterType.STRING,
                    description="Output file format",
                    required=False,
                    default="png",
                    enum=["png", "svg", "html"],
                ),
            ],
            category="drugdiscovery",
        )

    def check_available(self) -> bool:
        try:
            from rdkit import Chem  # noqa: F401
            from rdkit.Chem import Draw  # noqa: F401

            return True
        except ImportError:
            return False

    def _render_2d_png(self, mol, output_path: Path) -> None:
        """Render 2D structure as PNG."""
        from rdkit.Chem import Draw

        img = Draw.MolToImage(mol, size=(600, 400))
        img.save(str(output_path))

    def _render_2d_svg(self, mol, output_path: Path) -> None:
        """Render 2D structure as SVG."""
        from rdkit.Chem.Draw import rdMolDraw2D

        drawer = rdMolDraw2D.MolDraw2DSVG(600, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        output_path.write_text(svg)

    def _render_3d_html(self, mol, smiles: str, output_path: Path) -> None:
        """Generate an interactive 3D HTML viewer."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Generate 3D coordinates
        mol_3d = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol_3d, AllChem.ETKDGv3())
        if result != 0:
            # Fallback: use random coordinates
            AllChem.EmbedMolecule(mol_3d, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol_3d, maxIters=500)

        mol_block = Chem.MolToMolBlock(mol_3d)

        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>3D Molecule Viewer - {smiles}</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; padding: 20px; font-family: Arial, sans-serif; }}
        #viewer {{ width: 800px; height: 600px; position: relative; border: 1px solid #ccc; }}
        .info {{ margin-top: 10px; color: #555; }}
    </style>
</head>
<body>
    <h3>3D Molecule Viewer</h3>
    <p class="info">SMILES: <code>{smiles}</code></p>
    <div id="viewer"></div>
    <script>
        var viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
        var molData = `{mol_block}`;
        viewer.addModel(molData, "sdf");
        viewer.setStyle({{}}, {{stick: {{radius: 0.15}}, sphere: {{scale: 0.25}}}});
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>"""
        output_path.write_text(html_content)

    async def execute(self, **kwargs) -> ToolResult:
        smiles: str = kwargs["smiles"]
        mode: str = kwargs.get("mode", "2d")
        output_format: str = kwargs.get("output_format", "png")

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem  # noqa: F401
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

        # Compute 2D coordinates for 2D rendering
        if mode == "2d":
            AllChem.Compute2DCoords(mol)

        # Determine output path
        settings = get_settings()
        output_dir = settings.data_dir / "visualizations"
        output_dir.mkdir(parents=True, exist_ok=True)

        file_id = uuid.uuid4().hex[:8]

        # Validate format/mode combinations
        if mode == "3d" and output_format != "html":
            output_format = "html"

        ext = output_format
        output_path = output_dir / f"mol_{file_id}.{ext}"

        warnings: list[str] = []
        try:
            if mode == "2d":
                if output_format == "png":
                    self._render_2d_png(mol, output_path)
                elif output_format == "svg":
                    self._render_2d_svg(mol, output_path)
                elif output_format == "html":
                    # For 2D HTML, embed SVG
                    self._render_2d_svg(mol, output_path.with_suffix(".svg"))
                    output_path = output_path.with_suffix(".svg")
                    warnings.append("2D HTML not supported; generated SVG instead.")
            else:
                self._render_3d_html(mol, smiles, output_path)
        except Exception as e:
            return ToolResult(
                success=False,
                error=f"Rendering failed: {e}",
                summary=f"Failed to render molecule: {e}",
            )

        return ToolResult(
            success=True,
            data={
                "smiles": smiles,
                "mode": mode,
                "format": output_format,
                "output_path": str(output_path),
            },
            summary=f"Rendered {mode.upper()} visualization of molecule to {output_path}.",
            files=[output_path],
            warnings=warnings,
        )
