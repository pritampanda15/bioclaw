# BioClaw

Bioinformatics AI agent platform for drug discovery and NGS analysis. BioClaw combines a Claude-powered ReAct agent with a curated toolkit of 17 bioinformatics tools, accessible through a CLI, REST API, or Next.js web interface.

## Features

- **AI Agent** — Conversational interface backed by Claude that reasons through multi-step bioinformatics workflows
- **17 Built-in Tools** — Drug discovery (molecular properties, docking, QSAR, ADMET) and NGS (QC, alignment, variant calling, DE analysis)
- **Multiple Interfaces** — Interactive CLI, FastAPI server with WebSocket chat, and a Next.js 15 frontend with 3D molecule visualization
- **Extensible Tool System** — Add tools by dropping a Python file into a category directory; auto-discovered at startup
- **Background Jobs** — Long-running tools (docking, alignment, pipelines) execute asynchronously with status tracking

## Quick Start

### Installation

```bash
# Core install
pip install -e .

# With drug discovery tools (RDKit, Meeko, etc.)
pip install -e ".[drugdiscovery]"

# With NGS tools (pysam, PyDESeq2)
pip install -e ".[ngs]"

# Everything
pip install -e ".[all]"

# Development
pip install -e ".[dev]"
```

Some tools require system binaries (not pip-installable):

```bash
conda install -c conda-forge -c bioconda vina samtools bcftools fastqc hisat2 subread
```

### Set your API key

```bash
export ANTHROPIC_API_KEY="sk-ant-..."
```

### Usage

```bash
# Interactive chat with the agent
bioclaw chat

# Run a tool directly
bioclaw run molecular_properties -P smiles "CC(=O)Oc1ccccc1C(=O)O"

# List available tools and check dependencies
bioclaw tools list
bioclaw tools check

# Start the API server
bioclaw serve

# Manage background jobs
bioclaw jobs list
bioclaw jobs status <job_id>
bioclaw jobs cancel <job_id>
```

## Configuration

BioClaw reads configuration from `bioclaw.yaml` in the working directory, environment variables with the `BIOCLAW_` prefix, or built-in defaults.

```yaml
model: claude-sonnet-4-5-20250929
max_tool_rounds: 15
data_dir: ./bioclaw_data

server:
  host: 127.0.0.1
  port: 8000

tools:
  enabled_categories:
    - drugdiscovery
    - ngs
    - common
```

## Tools

### Drug Discovery (9 tools)

| Tool | Description |
|------|-------------|
| `molecular_properties` | Calculate MW, LogP, TPSA, Lipinski Ro5, QED, PAINS alerts from SMILES |
| `compound_search` | Search ChEMBL or PubChem by name, synonym, or ID |
| `similarity_search` | Tanimoto similarity with ECFP4 or MACCS fingerprints |
| `mol_visualizer` | 2D/3D molecule rendering to PNG, SVG, or interactive HTML |
| `admet_predictor` | Predict absorption, distribution, metabolism, excretion, and toxicity |
| `scaffold_analysis` | Murcko scaffold extraction or R-group decomposition |
| `qsar_model` | Build Random Forest or XGBoost QSAR models from SMILES + activity data |
| `molecular_docking` | AutoDock Vina docking with automatic PDB-to-PDBQT receptor conversion |
| `pdb_to_pdbqt` | Clean PDB structures (remove waters/ions/ligands) and convert to PDBQT |

### NGS (6 tools)

| Tool | Description |
|------|-------------|
| `quality_control` | Run FastQC on FASTQ files |
| `bam_stats` | Compute BAM statistics (mapping rate, coverage, insert size) |
| `read_alignment` | Align reads with BWA-MEM or HISAT2, output sorted BAM |
| `variant_calling` | Call variants with bcftools or GATK |
| `diff_expression` | DESeq2 differential expression analysis via PyDESeq2 |
| `rnaseq_pipeline` | End-to-end RNA-seq: FastQC &rarr; HISAT2 &rarr; featureCounts |

### Common (2 tools)

| Tool | Description |
|------|-------------|
| `file_converter` | Convert between CSV, TSV, and JSON formats |
| `database_query` | Query the local BioClaw results database |

## Examples

### Molecular Docking (end-to-end from PDB)

```bash
bioclaw run molecular_docking \
  -P ligand_smiles "CC(=O)Oc1ccccc1C(=O)O" \
  -P receptor_file "protein.pdb" \
  -P center_x 15.0 -P center_y 20.0 -P center_z 25.0
```

The tool automatically cleans the PDB (removes waters, ions, ligands), converts to PDBQT, and runs Vina.

### ADMET Profiling via Chat

```
$ bioclaw chat
You: What are the ADMET properties of ibuprofen?
Agent: I'll calculate the ADMET properties for ibuprofen (CC(C)Cc1ccc(cc1)C(C)C(=O)O)...
```

### RNA-seq Pipeline

```bash
bioclaw run rnaseq_pipeline \
  -P reads_1 "sample_R1.fastq.gz" \
  -P reads_2 "sample_R2.fastq.gz" \
  -P reference "genome.fa" \
  -P annotation_gtf "genes.gtf" \
  -P output_dir "results/"
```

## Architecture

```
src/bioclaw/
├── agent/          # ReAct loop, Claude client, conversation memory
├── cli/            # Click CLI (chat, run, serve, tools, jobs)
├── config/         # Pydantic Settings + YAML overlay
├── server/         # FastAPI app, WebSocket, REST routes
├── skills/         # Tool system
│   ├── base.py     # BaseTool ABC
│   ├── types.py    # ToolMetadata, ToolParameter, ToolResult
│   ├── registry.py # Auto-discovery via pkgutil
│   ├── drugdiscovery/
│   ├── ngs/
│   └── common/
├── storage/        # SQLite-backed results & job persistence
└── utils/          # Logging, validators

web/                # Next.js 15 frontend (React 19, Tailwind 4)
├── src/app/        # Pages: chat, jobs, settings
└── src/components/ # Chat, molecule viewer (3Dmol.js), NGS viz (Plotly)
```

**Agent loop** — The ReAct agent iterates up to 15 rounds of Think &rarr; Tool Use &rarr; Observe, streaming events (text deltas, tool calls, results) to the UI.

**Tool system** — Each tool is a `BaseTool` subclass with a `metadata` property and async `execute()`. Tools are auto-discovered by scanning category directories at startup. The registry checks binary and library availability before execution.

**Web frontend** — Next.js proxies `/api` to the FastAPI backend. The chat interface uses WebSockets for real-time streaming. Molecule results render with 3Dmol.js; NGS results with Plotly.

## Adding a Custom Tool

Create a file in the appropriate category directory (e.g., `src/bioclaw/skills/drugdiscovery/my_tool.py`):

```python
from bioclaw.skills.base import BaseTool
from bioclaw.skills.types import ParameterType, ToolMetadata, ToolParameter, ToolResult

class MyTool(BaseTool):
    @property
    def metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="my_tool",
            description="Description for the agent",
            parameters=[
                ToolParameter(
                    name="input_data",
                    type=ParameterType.STRING,
                    description="The input to process",
                ),
            ],
            category="drugdiscovery",
        )

    async def execute(self, **kwargs) -> ToolResult:
        data = kwargs["input_data"]
        # ... your logic ...
        return ToolResult(success=True, data={"result": data}, summary="Done")
```

The tool will be discovered automatically on next startup.

## Web Frontend

```bash
cd web
npm install
npm run dev
```

Open `http://localhost:3000`. Make sure the API server is running (`bioclaw serve`) on port 8000.

## Development

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

