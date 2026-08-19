# SiGMA: Single-cell Genome sequencing: Mechanistic Amplification simulator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**SiGMA** is a computational framework for simulating single-cell DNA sequencing (scDNA-seq) data. It first models somatic evolution along a coalescent cell lineage tree to generate ground-truth single-cell genomes alongside matched bulk profiles (including normal, tumor, and lineage-specific kindred samples). Next, SiGMA performs a mechanistic, cycle-by-cycle simulation of MDA and PTA by explicitly modeling hexamer priming, amplicon extension, and the accumulation of amplification bias across cycles to produce protocol-specific sequencing reads. This allows users to generate highly realistic, customizable synthetic datasets with known ground truths for SNVs.

This repository contains the source code, installation instructions, and tutorials required to work with SiGMA.

---

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Quick Start / Demo](#quick-start--demo)
- [Detailed Usage](#detailed-usage)
- [Output Description](#output-description)
- [Citation](#citation)
- [Contact](#contact)

---

## System Requirements

### Hardware requirements
SiGMA requires a standard computer with enough RAM to support the in-memory operations of genome processing. For optimal performance, we recommend:
*   **RAM:** 16+ GB
*   **CPU:** 4+ cores
*   **Storage:** 50+ GB available space (depends heavily on the size of the genome and the number of simulated cells)

### Software requirements
SiGMA is supported for *Linux* and *macOS*. The software has been tested on the following system:
*   22.04.5 LTS (GNU/Linux 6.8.0-111-generic x86_64)

**Dependencies:**
*   Python >= 3.8
---

## Installation
**1. Clone the repository**
First, download the SiGMA repository to your local machine to access the environment files and scripts:
```bash
git clone https://github.com/Sjd-Bh/SiGMA.git
cd SiGMA
```

**2. Create the Conda environments**
SiGMA relies on two specific environments to prevent package conflicts. Create them using the provided `.yml` files:
```bash
conda env create -f SingleCellSim_env.yml
conda env create -f picard_env.yml
`
## Usage

SiGMA is run through a single wrapper script, `sigma.sh`, which automatically handles environment switching and path resolution.

The pipeline sequentially executes three main simulation steps:
1. **Step 1: Coalescent Tree Making** 
2. **Step 2: Single-cell DNA Amplification and Sequencing Simulation** 
3. **Step 3: Bulk Simulation** 


You can install SiGMA by cloning this repository and installing the required dependencies. We recommend using a virtual environment (like `conda` or `venv`).
```bash
# Clone the repository
git clone https://github.com/Sjd-Bh/SiGMA.git
cd SiGMA

# Create and activate a conda environment
conda create -n sigma_env python=3.9
conda activate sigma_env

# Install dependencies
pip install -r requirements.txt
*Typical installation time on a standard desktop computer: 2-5 minutes.*

---

## Quick Start / Demo

To verify that SiGMA is installed correctly, you can run the included demo dataset. This demo simulates a small population of 10 cells based on a subset of the human genome (e.g., chr22).

1. Navigate to the demo directory:
bash
cd demo

2. Run the simulation script:
bash
python ../sigma.py --config demo_config.yaml --outdir ./demo_output

**Expected Output:**
The script should finish without errors and generate the following files in the `demo_output` folder:
*   `simulated_cells.bam`
*   `ground_truth_cnv.tsv`
*   `run_summary.log`

*Expected run time for demo on a standard desktop computer: ~10 minutes.*

---

## Detailed Usage

To run SiGMA on your own data, you need to provide a reference genome and a configuration file specifying the parameters of the simulation (e.g., number of cells, mutation rate, sequencing depth).

bash
python sigma.py -r [REFERENCE.fa] -c [CONFIG.yaml] -o [OUTPUT_DIR]

### Key Parameters:
*   `-r`, `--reference`: Path to the reference genome FASTA file.
*   `-c`, `--config`: Path to the YAML configuration file containing simulation parameters.
*   `-n`, `--num_cells`: Number of single cells to simulate (Overrides config file).
*   `-d`, `--depth`: Mean sequencing depth per cell.
*   `-t`, `--threads`: Number of CPU threads to use.

Please refer to the `docs/` folder for a comprehensive guide on configuring the YAML file.

---

## Output Description

SiGMA generates several output files designed to integrate seamlessly into existing bioinformatics pipelines:

1.  **`[PREFIX]_R1.fastq.gz` / `[PREFIX]_R2.fastq.gz`**: Simulated paired-end sequencing reads.
2.  **`ground_truth_events.bed`**: A BED file containing the exact genomic coordinates of all simulated SVs and CNVs.
3.  **`cell_lineage_tree.newick`**: The phylogenetic tree representing the clonal evolution of the simulated cells.

---

## Citation

If you use SiGMA in your research, please cite our paper:

> **[First Author Last Name] et al.** "SiGMA: A highly realistic single-cell DNA sequencing simulator." *Cell Reports Methods* (202X). DOI: [Insert DOI here once published/accepted]

---

## Contact

For bug reports, feature requests, or general questions, please open an issue on the [GitHub Issues](https://github.com/Sjd-Bh/SiGMA/issues) page or contact [Your Name/Email].

