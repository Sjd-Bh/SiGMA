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
```

## Usage

SiGMA is run through a single wrapper script, `sigma.sh`, which automatically handles environment switching and path resolution.

The pipeline sequentially executes three main simulation steps:
1. **Step 1: Coalescent Tree Making** 
2. **Step 2: Single-cell DNA Amplification and Sequencing Simulation** 
3. **Step 3: Bulk Simulation** 


## Quick Demo

To verify that SiGMA is installed correctly, you can run a quick simulation using the 1Mb test reference file provided in the repository. 

Make sure you are in the main `SiGMA` directory, then run the following command:
```bash
# Make the script executable (only needed once)
chmod +x sigma.sh
```

### Step 1: Coalescent Tree Simulation
This step generates the evolutionary tree for the cells. 
```bash
./sigma.sh coalescent \
  --outdir test_output \
  --n_cells 10 \
  --genome_length 1000000
```

### Step 2: Single-Cell Sequencing Simulation
This step generates the single-cell genomes, introduces mutations, simulates amplification (MDA/PTA), and generates downsampled BAM files. *(Replace `/path/to/picard.jar` with your actual path).*
```bash
./sigma.sh scDNAseq \
  --ref test/ref.fa \
  --picard_jar /path/to/picard.jar \
  --outdir test_output
```

### Step 3: Bulk DNA Sequencing Simulation
This final step simulates matched normal and tumor bulk sequencing FASTQ and BAM files based on the cell outputs from the previous steps.
```bash
./sigma.sh bulk \
  --ref test/ref.fa \
  --outdir test_output
```

**Expected Output:**
Once all three steps are completed, your `test_output/` directory will be populated with:
* `coal/`: The simulated coalescent tree (`coalTree.pkl`).
* `sc/`: Single-cell mutated genomes, variant VCFs, amplification directories, and final single-cell BAM files.
* `bulk/`: Simulated matched normal and tumor bulk FASTQ and BAM files.

---

## Detailed Usage

The `sigma.sh` wrapper script is the main entry point for the pipeline. It is executed using the following syntax:
```bash
./sigma.sh {coalescent|scDNAseq|bulk} [OPTIONS]
```

### Commands
SiGMA is divided into three sequential subcommands:
*   `coalescent`: Simulates the evolutionary coalescent tree of the cells.
*   `scDNAseq`: Simulates single-cell genomes, introduces SNVs and CNVs, models amplification (MDA/PTA), and generates downsampled single-cell BAM files.
*   `bulk`: Simulates matched normal and tumor bulk sequencing FASTQ and BAM files.

---

### Essential Arguments
These arguments **must** be provided depending on the specific subcommand you are running.

*   `--ref PATH`
*   **Description:** Path to the reference FASTA file.
*   **Required for:** `scDNAseq` and `bulk` commands.
*   `--picard_jar PATH`
*   **Description:** Path to the downloaded `picard.jar` executable.
*   **Required for:** `scDNAseq` command.

---

### Optional Arguments
These arguments have pre-configured defaults, allowing you to run SiGMA out-of-the-box. You can override them to customize your simulation.

#### General Options
*   `--outdir PATH`
*   **Description:** Directory where all simulation outputs will be saved.
*   **Default:** `./sigma_output` (in your current working directory)
*   `--cores` or `--num_cores INT`
*   **Description:** Number of CPU cores to allocate for parallel processing steps.
*   **Default:** `10`

#### Tree & Genome Parameters (Used in `coalescent` & `scDNAseq`)
*   `--n_cells INT`
*   **Description:** Total number of cells to simulate in the coalescent tree.
*   **Default:** `10`
*   `--genome_length INT`
*   **Description:** The length of the genomic region to simulate (in base pairs).
*   **Default:** `1000000` (1Mb)
*   `--chrom STR`
*   **Description:** Chromosome name used in the generated VCFs and FASTAs.
*   **Default:** Dynamically calculated based on genome length (e.g., `1000kb` for a 1,000,000 bp genome length).

#### Mutation Parameters (Used in `scDNAseq`)
*   `--mutation_rate FLOAT`
*   **Description:** Single Nucleotide Variant (SNV) somatic mutation rate per base pair.
*   **Default:** `1e-6`
*   `--cnv_rate FLOAT`
*   **Description:** Copy Number Variation (CNV) rate.
*   **Default:** `4e-9`
*   `--mean_cnv_length INT`
*   **Description:** Mean length of generated Copy Number Variations.
*   **Default:** `1000`

#### Single-Cell & Amplification Parameters (Used in `scDNAseq`)
*   `--node INT`
*   **Description:** Specific tree node (cell) to extract for downstream single-cell amplification simulation.
*   **Default:** A random node chosen between `0` and `n_cells - 1`.
*   `--num_simulations INT`
*   **Description:** Number of independent amplification (MDA/PTA) simulations to run.
*   **Default:** `10`
*   `--target_depth INT`
*   **Description:** Target sequencing depth for BQSR downsampling.
*   **Default:** Automatically calculated based on the minimum mean depth found in intermediate BAMs (falls back to `26` if auto-detection fails).

#### Bulk Sequencing Parameters (Used in `bulk`)
*   `--tumor_coverage FLOAT`
*   **Description:** Target sequencing coverage for the simulated tumor bulk BAM.
*   **Default:** `40.0`
*   `--normal_coverage FLOAT`
*   **Description:** Target sequencing coverage for the simulated normal bulk BAM.
*   **Default:** `40.0`
*   `--read_length INT`
*   **Description:** Read length for the generated paired-end bulk FASTQ files.
*   **Default:** `150`
*   `--insert_size INT`
*   **Description:** Mean insert size for the bulk paired-end reads.
*   **Default:** `300`
*   `--std_dev INT`
*   **Description:** Standard deviation for the bulk read insert size.
*   **Default:** `20`

---

## Citation



---

## Contact

For bug reports, feature requests, or general questions, please open an issue on the [GitHub Issues](https://github.com/Sjd-Bh/SiGMA/issues) page or contact [sajedeh.bahonar@gmail.com].

