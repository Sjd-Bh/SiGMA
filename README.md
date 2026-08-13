# SiGMA: [Single-cell Genome sequencing: Mechanistic Amplification simulator<img width="468" height="13" alt="image" src="https://github.com/user-attachments/assets/3f053f1f-64d7-4257-809f-5c30cf9fd7fd" />
]

**SiGMA** is a computational tool for simulating single-cell DNA sequencing (scDNA-seq) data. It allows users to generate highly realistic, customizable synthetic datasets with known ground truths for [copy number variations, structural variants, read depth overdispersion, etc.], enabling the robust benchmarking of downstream bioinformatics tools.

This repository contains the source code, installation instructions, and tutorials required to reproduce the results presented in our *Cell Reports Methods* publication.

---

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Quick Start / Demo](#quick-start--demo)
- [Detailed Usage](#detailed-usage)
- [Output Description](#output-description)
- [Reproducing Paper Results](#reproducing-paper-results)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## System Requirements

### Hardware requirements
SiGMA requires a standard computer with enough RAM to support the in-memory operations of genome processing. For optimal performance, we recommend:
*   **RAM:** 16+ GB
*   **CPU:** 4+ cores
*   **Storage:** 50+ GB available space (depends heavily on the size of the genome and the number of simulated cells)

### Software requirements
SiGMA is supported for *Linux* and *macOS*. The software has been tested on the following systems:
*   Ubuntu 20.04 / 22.04
*   macOS Monterey (12.0)

**Dependencies:**
*   Python >= 3.8
*   [List other dependencies, e.g., samtools >= 1.10, bedtools >= 2.29]
*   [List key Python packages: numpy, pandas, pysam, scipy]

---

## Installation

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

## Reproducing Paper Results

To reproduce the exact figures and benchmarking results presented in our *Cell Reports Methods* paper, please navigate to the `paper_analysis/` directory. 

bash
cd paper_analysis
bash run_all_benchmarks.sh
*Note: Reproducing the full paper results requires downloading the complete hg38 reference genome and may take several days on a high-performance compute cluster.* Detailed instructions are available in `paper_analysis/README.md`.

---

## Citation

If you use SiGMA in your research, please cite our paper:

> **[First Author Last Name] et al.** "SiGMA: A highly realistic single-cell DNA sequencing simulator." *Cell Reports Methods* (202X). DOI: [Insert DOI here once published/accepted]

---

## License

This project is covered under the **MIT License**. See the `LICENSE` file for more details.

---

## Contact

For bug reports, feature requests, or general questions, please open an issue on the [GitHub Issues](https://github.com/Sjd-Bh/SiGMA/issues) page or contact [Your Name/Email].


### How to use this:
1. Create a file named `README.md` on your server in the root of your `SiGMA` folder (`nano README.md`).
2. Paste this text in.
3. Edit the bracketed placeholders to match your exact tool's parameters.
4. Save the file.
5. Push the new README to GitHub:
