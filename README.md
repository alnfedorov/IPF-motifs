This repository contains the code and data to reproduce the motif analysis presented in the paper "....".

## Repository Structure

```
IPF-motifs/
├── GRCh38/     # GRCh38 reference genome files and processed annotations (populated by setup)
├── pixi.lock   # Pixi lock file for reproducible environments
├── pixi.toml   # Pixi project configuration, dependencies, and tasks
├── setup/      # Scripts for initial data download and setup
├── utils/      # Utility scripts and modules (e.g., for FASTA, BED, motif operations)
└── stories/    # Main analysis workflows
    ├── cCRE/   # ENCODE cCRE processing and matching to transcripts
    ├── motifs/ # Motif database parsing, promoter scoring, and per-gene response calculation
    └── IPF/    # IPF scRNA-seq data integration and significance testing
```

## Installation & Setup

1. **Install Pixi:** The computational environment and all dependencies are managed by Pixi. Install it
   from [pixi.sh](https://pixi.sh/).
2. **Clone the Repository:**
   ```bash
   git clone https://github.com/alnfedorov/IPF-motifs
   cd IPF-motifs
   ```
3. **Download and Prepare Genomic Data:**
   ```bash
   pixi run setup
   ```
   This task downloads the GRCh38 reference genome and GENCODE annotations, then processes them. ENCODE cCREs, JASPAR
   motifs, and scRNA-seq tables are pre-committed to this repository within the relevant `GRCh38/` or
   `stories/*/ld/resources/` directories.

**Note:** This repository is designed to be used on Linux x86_64 systems and is not compatible with other architectures
or operating systems.

## Workflow

The analysis is divided into stages called 'stories'. Each story corresponds to a Pixi task and stores its results in a
local `ld/results/` directory. The `ld/` (local data) subdirectories also contain story-specific configurations (e.g.,
file paths, thresholds) and input resources.

Execute the following Pixi tasks in order to reproduce the full analysis:

1. **Process cCREs:**
   This stage identifies promoter-like sequences (PLS) from ENCODE cCREs, matches them to protein-coding gene
   transcripts, imputes promoters for transcripts lacking ENCODE PLS, normalizes their lengths, and extracts their DNA
   sequences.
   ```bash
   pixi run stories/cCRE
   ```
2. **Score Motifs & Calculate Per-Gene Responses:**
   Promoter sequences are scanned using non-redundant JASPAR motifs. Scores are then aggregated to calculate a response
   value for each gene-motif cluster combination.
   ```bash
   pixi run stories/motifs
   ```
3. **IPF-Specific Analysis & Plotting:**
   Motif cluster responses are integrated with IPF scRNA-seq expression data to identify statistically significant
   associations.
   ```bash
   pixi run stories/IPF
   ```

## Dependencies

All Python and system tool dependencies are defined in `pixi.toml` and managed by Pixi. They will be automatically
installed into an isolated environment when you run any Pixi command (like `pixi run <task>`) or explicitly with
`pixi install`.
