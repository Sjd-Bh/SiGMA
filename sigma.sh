#!/bin/bash
# sigma.sh - Master Wrapper for SiGMA Simulator

set -euo pipefail

# ---------------------------------------------------------
# Resolve SiGMA Directory (Allows running from anywhere)
# ---------------------------------------------------------
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SIGMA_DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

# Let Python find SinglCellSim when runAmp.py is launched from this wrapper
export PYTHONPATH="$SIGMA_DIR${PYTHONPATH:+:$PYTHONPATH}"

# ---------------------------------------------------------
# Help Menu Function
# ---------------------------------------------------------
print_help() {
    echo "Usage: bash sigma.sh {coalescent|scDNAseq|bulk} [OPTIONS]"
    echo ""
    echo "Commands:"
    echo "  coalescent    Run Coalescent Tree Simulation"
    echo "  scDNAseq      Run Single Cell Sequencing Simulation"
    echo "  bulk          Run Bulk DNAseq Simulation"
    echo ""
    echo "Options:"
    echo "  -h, --help            Show this help message and exit"
    echo "  --outdir PATH         Directory to save outputs (Default: ./sigma_output)"
    echo "  --ref PATH            Path to reference fasta (Required for scDNAseq and bulk)"
    echo "  --picard_jar PATH     Path to picard.jar (Required for scDNAseq)"
    echo "  --amp {MDA|PTA}       Amplification method (Default: MDA)"
    echo "                        MDA -> SinglCellSim/configs/config_files/MDAsim.ini"
    echo "                        PTA -> SinglCellSim/configs/config_files/PTAsim.ini"
    echo "  --n_cells INT         Number of cells (Default: 10)"
    echo "  --genome_length INT   Length of genome (Default: 1000000)"
    echo "  --chrom STR           Chromosome name (Default: calculated based on length, e.g. 1000kb)"
    echo "  --mutation_rate FLOAT Mutation rate (Default: 1e-6)"
    echo "  --cnv_rate FLOAT      CNV rate (Default: 4e-9)"
    echo "  --mean_cnv_length INT Mean CNV length (Default: 1000)"
    echo "  --node INT            Tree node to extract (Default: Random node)"
    echo "  --num_simulations INT Number of simulations (Default: 10)"
    echo "  --cores, --num_cores  Number of CPU cores (Default: 10)"
    echo "  --target_depth INT    Target depth for downsampling (Default: Auto-extracted minimum)"
    echo "  --tumor_coverage FLT  Tumor coverage for bulk (Default: 40.0)"
    echo "  --normal_coverage FLT Normal coverage for bulk (Default: 40.0)"
    echo "  --read_length INT     Read length (Default: 150)"
    echo "  --insert_size INT     Insert size (Default: 300)"
    echo "  --std_dev INT         Standard deviation (Default: 20)"
    echo ""
}

# Check for help flag or empty arguments
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || -z "${1:-}" ]]; then
    print_help
    exit 0
fi

COMMAND=$1
shift

# ---------------------------------------------------------
# Default values
# ---------------------------------------------------------
OUT_DIR="$PWD/sigma_output"
REF=""
PICARD_JAR=""
AMP="MDA"
N_CELLS=10
GENOME_LENGTH=1000000
CHROM=""
MUT_RATE="1e-6"
CNV_RATE="4e-9"
MEAN_CNV_LEN="1000"
NODE=""
NUM_SIMS="10"
CORES="10"
TARGET_DEPTH=""
TUMOR_COV="40.0"
NORMAL_COV="40.0"
READ_LEN="150"
INSERT_SIZE="300"
STD_DEV="20"

# ---------------------------------------------------------
# Parse remaining arguments
# ---------------------------------------------------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --outdir) OUT_DIR="$2"; shift ;;
        --ref) REF="$2"; shift ;;
        --picard_jar) PICARD_JAR="$2"; shift ;;
        --amp) AMP="$2"; shift ;;
        --n_cells) N_CELLS="$2"; shift ;;
        --genome_length) GENOME_LENGTH="$2"; shift ;;
        --chrom) CHROM="$2"; shift ;;
        --mutation_rate) MUT_RATE="$2"; shift ;;
        --cnv_rate) CNV_RATE="$2"; shift ;;
        --mean_cnv_length) MEAN_CNV_LEN="$2"; shift ;;
        --node) NODE="$2"; shift ;;
        --num_simulations) NUM_SIMS="$2"; shift ;;
        --num_cores|--cores) CORES="$2"; shift ;;
        --target_depth) TARGET_DEPTH="$2"; shift ;;
        --tumor_coverage) TUMOR_COV="$2"; shift ;;
        --normal_coverage) NORMAL_COV="$2"; shift ;;
        --read_length) READ_LEN="$2"; shift ;;
        --insert_size) INSERT_SIZE="$2"; shift ;;
        --std_dev) STD_DEV="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Normalize and validate amplification method
AMP="$(echo "$AMP" | tr '[:lower:]' '[:upper:]')"
case "$AMP" in
    MDA|PTA) ;;
    *)
        echo "Error: --amp must be MDA or PTA (got: $AMP)"
        exit 1
        ;;
esac

# ---------------------------------------------------------
# Dynamic Defaults Logic
# ---------------------------------------------------------
if [ -z "$CHROM" ]; then
    CHROM="$((${GENOME_LENGTH} / 1000))kb"
fi

if [ -z "$NODE" ]; then
    NODE=$((RANDOM % N_CELLS))
fi

# ---------------------------------------------------------
# Define and create base output directories
# ---------------------------------------------------------
COAL_OUT="$OUT_DIR/coal"
SC_OUT="$OUT_DIR/sc"
BULK_OUT="$OUT_DIR/bulk"
AMP_OUT="$SC_OUT/amp/$AMP"

mkdir -p "$COAL_OUT"
mkdir -p "$SC_OUT/snps"
mkdir -p "$SC_OUT/mutations/cells_output"
mkdir -p "$SC_OUT/cell_genome"
mkdir -p "$SC_OUT/amp/MDA"
mkdir -p "$SC_OUT/amp/PTA"
mkdir -p "$SC_OUT/bam"
mkdir -p "$BULK_OUT/fastqFiles"
mkdir -p "$BULK_OUT/bamfiles"

resolve_amp_config() {
    local candidates=(
        "$SIGMA_DIR/SinglCellSim/configs/config_files/${AMP}sim.ini"
        "$SIGMA_DIR/SinglCellSim/configs/${AMP}sim.ini"
    )
    local cfg
    for cfg in "${candidates[@]}"; do
        if [ -f "$cfg" ]; then
            echo "$cfg"
            return 0
        fi
    done
    echo "Error: ${AMP} config not found. Looked for:" >&2
    printf '  %s\n' "${candidates[@]}" >&2
    return 1
}

# Fallback if real_data/finding_mean_depth.sh is absent from the checkout
compute_mean_depths() {
    local bam_list="$1"
    local out_csv="$2"
    echo "bam,mean_depth" > "$out_csv"
    while IFS= read -r bam || [ -n "$bam" ]; do
        [ -z "$bam" ] && continue
        if [ ! -f "$bam" ]; then
            echo "Warning: BAM not found, skipping: $bam" >&2
            continue
        fi
        local depth
        depth="$(samtools depth "$bam" | awk '{s+=$3; n++} END{if(n) printf "%.4f", s/n; else print 0}')"
        echo "${bam},${depth}" >> "$out_csv"
    done < "$bam_list"
}

# ---------------------------------------------------------
# Main Execution Blocks
# ---------------------------------------------------------
case "$COMMAND" in

  coalescent)
    echo "=== Running Coalescent Tree Simulation ==="
    python "$SIGMA_DIR/SinglCellSim/Coal/coalescentSim.py" \
      --n_cells "$N_CELLS" --genome_length "$GENOME_LENGTH" --N 100 \
      --output "$COAL_OUT/coalTree.pkl"

    echo "Coalescent simulation complete. Outputs saved to $COAL_OUT"
    ;;

  scDNAseq)
    if [ -z "$REF" ] || [ -z "$PICARD_JAR" ]; then
      echo "Error: Both --ref and --picard_jar are required for scDNAseq."
      exit 1
    fi

    AMP_CONFIG="$(resolve_amp_config)"
    BAM_LIST="$AMP_OUT/${AMP}_bam.txt"
    DEPTH_CSV="$SC_OUT/bam/mean_depth_${AMP}_bam.csv"
    DOWNSAMPLED_LIST="$AMP_OUT/downsampled_${AMP}_bam.txt"

    echo "=== Running Single Cell Sequencing Simulation ==="
    echo "[INFO] Amplification method: $AMP"
    echo "[INFO] Amp config: $AMP_CONFIG"
    echo "[INFO] Amp output: $AMP_OUT"

    # 2. SNPs
    python3 "$SIGMA_DIR/SinglCellSim/generateCellGenme/simSNPs.py" \
      -o "$SC_OUT/snps" \
      -R "$REF" \
      --rate 0.001

    python3 "$SIGMA_DIR/SinglCellSim/generateCellGenme/introduceSNPs.py" \
      --ref "$REF" \
      --vcf "$SC_OUT/snps/paternal_snps.vcf" \
      --out "$SC_OUT/cell_genome/paternal_genome.fasta" \
      --chrom "$CHROM"

    python3 "$SIGMA_DIR/SinglCellSim/generateCellGenme/introduceSNPs.py" \
      --ref "$REF" \
      --vcf "$SC_OUT/snps/maternal_snps.vcf" \
      --out "$SC_OUT/cell_genome/maternal_genome.fasta" \
      --chrom "$CHROM"

    # 3 & 4. Tree mutations and Cell Genome
    python "$SIGMA_DIR/SinglCellSim/Coal/treeMutations.py" \
      --input "$COAL_OUT/coalTree.pkl" \
      --output "$SC_OUT/mutations/mutations_coalTree.pkl" \
      --genome_length "$GENOME_LENGTH" \
      --mutation_rate "$MUT_RATE" \
      --cnv_rate "$CNV_RATE" \
      --mean_cnv_length "$MEAN_CNV_LEN" \
      --cell_out_dir "$SC_OUT/mutations/cells_output/" \
      --patref "$SC_OUT/cell_genome/paternal_genome.fasta" \
      --matref "$SC_OUT/cell_genome/maternal_genome.fasta" \
      --chrom "$CHROM" \
      --pat_prob 0.5

    python "$SIGMA_DIR/SinglCellSim/generateCellGenme/applySNVsCNVs_toCellGenome.py" \
      --pat "$SC_OUT/cell_genome/paternal_genome.fasta" \
      --mat "$SC_OUT/cell_genome/maternal_genome.fasta" \
      --pkl "$SC_OUT/mutations/mutations_coalTree.pkl" \
      --node "$NODE" \
      --output "$SC_OUT/cell_genome/"

    # Use the node-specific CNV BED written by treeMutations.py
    CNV_BED="$SC_OUT/mutations/cells_output/${NODE}_cnvs.bed"
    if [ ! -f "$CNV_BED" ]; then
      echo "Error: CNV BED not found: $CNV_BED"
      ls -l "$SC_OUT/mutations/cells_output/"
      exit 1
    fi
    echo "[INFO] Using CNV BED: $CNV_BED"

    # 5. Amplification: MDA -> MDAsim.ini, PTA -> PTAsim.ini
    python "$SIGMA_DIR/SinglCellSim/AmpSim/runAmp.py" \
      --num_simulations "$NUM_SIMS" \
      --config_file "$AMP_CONFIG" \
      --output_base "$AMP_OUT/" \
      --patSeq_file "$SC_OUT/cell_genome/${NODE}_paternal.fasta" \
      --matSeq_file "$SC_OUT/cell_genome/${NODE}_maternal.fasta" \
      --cnv_bed_file "$CNV_BED" \
      --num_cores "$CORES"

    shopt -s nullglob
    subset_fastas=("$AMP_OUT"/sim*/subset.fasta)
    shopt -u nullglob
    if [ ${#subset_fastas[@]} -eq 0 ]; then
      echo "Error: no subset.fasta files found under $AMP_OUT/sim*/"
      ls -ld "$AMP_OUT"/sim*/ || true
      exit 1
    fi
    echo "[INFO] Found ${#subset_fastas[@]} subset.fasta files for sequencing"

    # 6. Sequencing. Keep the glob as a string; scDNAseqSim.py expands it.
    python "$SIGMA_DIR/SinglCellSim/seqSim/scDNAseqSim.py" \
      --ref "$REF" \
      --scFiles "$AMP_OUT/sim*/subset.fasta" \
      --cores "$CORES"

    # Unquoted sim*/ glob is avoided by searching under AMP_OUT
    find "$AMP_OUT" -type f -name "*_subset.bam" | sort > "$BAM_LIST"
    if [ ! -s "$BAM_LIST" ]; then
      echo "Error: no *_subset.bam files found under $AMP_OUT"
      exit 1
    fi
    echo "[INFO] BAM list: $BAM_LIST"

    if [ -f "$SIGMA_DIR/real_data/finding_mean_depth.sh" ]; then
      bash "$SIGMA_DIR/real_data/finding_mean_depth.sh" "$BAM_LIST" -o "$SC_OUT/bam/"
      if [ ! -f "$DEPTH_CSV" ]; then
        if [ -f "$SC_OUT/bam/mean_depth_bam.csv" ]; then
          DEPTH_CSV="$SC_OUT/bam/mean_depth_bam.csv"
        else
          echo "Error: mean-depth CSV not written to $SC_OUT/bam/"
          ls -l "$SC_OUT/bam/" || true
          exit 1
        fi
      fi
    else
      echo "[WARN] real_data/finding_mean_depth.sh is missing; computing mean depths with samtools"
      compute_mean_depths "$BAM_LIST" "$DEPTH_CSV"
    fi

    if [ -z "$TARGET_DEPTH" ]; then
        TARGET_DEPTH="$(tail -n +2 "$DEPTH_CSV" | awk -F',' '{print $2}' | sort -n | head -n 1)"
        if [ -z "$TARGET_DEPTH" ]; then TARGET_DEPTH=26; fi
        echo "Auto-detected minimum target depth: $TARGET_DEPTH"
    fi

    bash "$SIGMA_DIR/real_data/downsam_BQ.sh" \
      --input "$BAM_LIST" \
      --picard "$PICARD_JAR" \
      --ref "$REF" \
      --dbsnp "$SC_OUT/snps/merged_snps.dedup.vcf.gz" \
      --csv "$DEPTH_CSV" \
      --target-depth "$TARGET_DEPTH"

    find "$AMP_OUT" -type f -name "*recal.bam" | sort > "$DOWNSAMPLED_LIST"
    if [ ! -s "$DOWNSAMPLED_LIST" ]; then
      echo "Error: no *recal.bam files found under $AMP_OUT"
      exit 1
    fi

    bash "$SIGMA_DIR/real_data/subset_sort_index.sh" \
      --input "$DOWNSAMPLED_LIST" \
      --cores "$CORES"

    echo "scDNAseq simulation complete. Outputs saved to $SC_OUT"
    ;;

  bulk)
    if [ -z "$REF" ]; then
      echo "Error: --ref is required for bulk"
      exit 1
    fi
    echo "=== Running Bulk DNAseq Simulation ==="

    python "$SIGMA_DIR/BulkSim/kindred_Normal_bulk_sim.py" \
      --paternal_fasta "$SC_OUT/cell_genome/paternal_genome.fasta" \
      --maternal_fasta "$SC_OUT/cell_genome/maternal_genome.fasta" \
      --cell_dir "$SC_OUT/mutations/cells_output/" \
      --output_dir "$BULK_OUT/fastqFiles/" \
      --tumor_coverage "$TUMOR_COV" \
      --normal_coverage "$NORMAL_COV" \
      --read_length "$READ_LEN" \
      --insert_size "$INSERT_SIZE" \
      --std_dev "$STD_DEV"

    python "$SIGMA_DIR/BulkSim/bulk_kindred_normal_sim_to_bam.py" \
      --ref "$REF" \
      --known_sites "$SC_OUT/snps/merged_snps.dedup.vcf.gz" \
      --r1 "$BULK_OUT/fastqFiles/normal_bulk_R1.fq" \
      --r2 "$BULK_OUT/fastqFiles/normal_bulk_R2.fq" \
      --sample_name Normal \
      --output_dir "$BULK_OUT/bamfiles/"

    python "$SIGMA_DIR/BulkSim/bulk_kindred_normal_sim_to_bam.py" \
      --ref "$REF" \
      --known_sites "$SC_OUT/snps/merged_snps.dedup.vcf.gz" \
      --r1 "$BULK_OUT/fastqFiles/match_tumor_bulk_R1.fq" \
      --r2 "$BULK_OUT/fastqFiles/match_tumor_bulk_R2.fq" \
      --sample_name match_tumor_bulk \
      --output_dir "$BULK_OUT/bamfiles/"

    parallel -j "$CORES" python "$SIGMA_DIR/BulkSim/bulk_kindred_normal_sim_to_bam.py" \
      --ref "$REF" \
      --known_sites "$SC_OUT/snps/merged_snps.sorted.vcf.gz" \
      --r1 "$BULK_OUT/fastqFiles/{}_clone_R1.fq" \
      --r2 "$BULK_OUT/fastqFiles/{}_clone_R2.fq" \
      --sample_name clone_{} \
      --output_dir "$BULK_OUT/bamfiles/" \
      --threads 8 ::: $(seq 0 $(($NUM_SIMS - 1)))

    echo "Bulk simulation complete. Outputs saved to $BULK_OUT"
    ;;

  *)
    echo "Invalid command: $COMMAND"
    echo "Run 'bash sigma.sh --help' for usage."
    exit 1
    ;;
esac
