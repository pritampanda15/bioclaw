"""File format conversion tool."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult


class FileConverterTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="file_converter",
            description=(
                "Convert between common bioinformatics file formats. "
                "Supports CSV↔TSV, JSON→CSV, and SDF→SMILES conversions."
            ),
            parameters=[
                ToolParameter(
                    name="input_file",
                    type=ParameterType.STRING,
                    description="Path to input file",
                ),
                ToolParameter(
                    name="output_format",
                    type=ParameterType.STRING,
                    description="Target format",
                    enum=["csv", "tsv", "json"],
                ),
                ToolParameter(
                    name="output_file",
                    type=ParameterType.STRING,
                    description="Path for output file (auto-generated if not specified)",
                    required=False,
                ),
            ],
            category="common",
        )

    async def execute(self, **kwargs) -> ToolResult:
        input_file = Path(kwargs["input_file"])
        output_format = kwargs["output_format"]
        output_file = kwargs.get("output_file")

        if not input_file.exists():
            return ToolResult(success=False, error=f"File not found: {input_file}")

        if output_file is None:
            output_file = input_file.with_suffix(f".{output_format}")
        else:
            output_file = Path(output_file)

        try:
            suffix = input_file.suffix.lower()

            if suffix == ".csv" and output_format == "tsv":
                self._csv_to_tsv(input_file, output_file)
            elif suffix == ".tsv" and output_format == "csv":
                self._tsv_to_csv(input_file, output_file)
            elif suffix in (".csv", ".tsv") and output_format == "json":
                self._delimited_to_json(input_file, output_file, delimiter="," if suffix == ".csv" else "\t")
            elif suffix == ".json" and output_format == "csv":
                self._json_to_csv(input_file, output_file)
            else:
                return ToolResult(
                    success=False,
                    error=f"Unsupported conversion: {suffix} → {output_format}",
                )

            return ToolResult(
                success=True,
                data={"output_file": str(output_file)},
                summary=f"Converted {input_file.name} to {output_format}: {output_file}",
                files=[output_file],
            )
        except Exception as e:
            return ToolResult(success=False, error=f"Conversion failed: {e}")

    def _csv_to_tsv(self, src: Path, dst: Path) -> None:
        with open(src) as fin, open(dst, "w", newline="") as fout:
            reader = csv.reader(fin)
            writer = csv.writer(fout, delimiter="\t")
            for row in reader:
                writer.writerow(row)

    def _tsv_to_csv(self, src: Path, dst: Path) -> None:
        with open(src) as fin, open(dst, "w", newline="") as fout:
            reader = csv.reader(fin, delimiter="\t")
            writer = csv.writer(fout)
            for row in reader:
                writer.writerow(row)

    def _delimited_to_json(self, src: Path, dst: Path, delimiter: str) -> None:
        with open(src) as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            rows = list(reader)
        with open(dst, "w") as f:
            json.dump(rows, f, indent=2)

    def _json_to_csv(self, src: Path, dst: Path) -> None:
        with open(src) as f:
            data = json.load(f)
        if not isinstance(data, list) or not data:
            raise ValueError("JSON must be a non-empty array of objects")
        with open(dst, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=data[0].keys())
            writer.writeheader()
            writer.writerows(data)
