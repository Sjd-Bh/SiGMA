#!/bin/bash
# sigma.sh - Master Wrapper for SiGMA Simulator

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
if [[ "$1" == "-h" || "$1" == "--help" || -z "$1" ]]; then
    print_help
    exit 0
fi

COMMAND=$1
shift # Shift arguments so $1 becomes the first flag

# ---------------------------------------------------------
# Default values
# ---------------------------------------------------------
OUT_DIR="$PWD/sigma_output"
REF=""
PICARD_JAR=""
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

mkdir -p "$COAL_OUT"
mkdir -p "$SC_OUT/snps"
mkdir -p "$SC_OUT/mutations/cells_output"
mkdir -p "$SC_OUT/cell_genome"
mkdir -p "$SC_OUT/amp/MDA"
mkdir -p "$SC_OUT/amp/PTA"
mkdir -p "$SC_OUT/bam"
mkdir -p "$BULK_OUT/fastqFiles"
mkdir -p "$BULK_OUT/bamfiles"

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
    echo "=== Running Single Cell Sequencing Simulation ==="
    
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
    
    # 5 & 6. Amp and scDNAseq Sim
    python "$SIGMA_DIR/SinglCellSim/AmpSim/runAmp.py" \
      --num_simulations "$NUM_SIMS" \
      --config_file "$SIGMA_DIR/SinglCellSim/configs/config_files/MDAsim.ini" \
      --output_base "$SC_OUT/amp/MDA/" \
      --patSeq_file "$SC_OUT/cell_genome/${NODE}_paternal.fasta" \
      --matSeq_file "$SC_OUT/cell_genome/${NODE}_maternal.fasta" \
      --cnv_bed_file "$SC_OUT/mutations/cnv.bed" \
      --num_cores "$CORES"
    
    python "$SIGMA_DIR/SinglCellSim/seqSim/scDNAseqSim.py" \
      --ref "$REF" \
      --scFiles "$SC_OUT/amp/PTA/sim*/subset.fasta" \
      --cores "$CORES"
    
    # Downsampling steps
    find "$SC_OUT/amp/MDA/sim*/" -type f -name "*_subset.bam" | sort > "$SC_OUT/amp/MDA/MDA_bam.txt"
    bash "$SIGMA_DIR/real_data/finding_mean_depth.sh" "$SC_OUT/amp/MDA/MDA_bam.txt" -o "$SC_OUT/bam/"
    
    # Determine Target Depth
    if [ -z "$TARGET_DEPTH" ]; then
        TARGET_DEPTH=$(tail -n +2 "$SC_OUT/bam/mean_depth_PTA_bam.csv" | awk -F',' '{print $2}' | sort -n | head -n 1)
        if [ -z "$TARGET_DEPTH" ]; then TARGET_DEPTH=26; fi
        echo "Auto-detected minimum target depth: $TARGET_DEPTH"
    fi

    # Conda environment handling
    #source $(conda info --base)/etc/profile.d/conda.sh
    #conda activate picard
    
    bash "$SIGMA_DIR/real_data/downsam_BQ.sh" \
      --input "$SC_OUT/amp/PTA/PTA_bam.txt" \
      --picard "$PICARD_JAR" \
      --ref "$REF" \
      --dbsnp "$SC_OUT/snps/merged_snps.dedup.vcf.gz" \
      --csv "$SC_OUT/bam/mean_depth_PTA_bam.csv" \
      --target-depth "$TARGET_DEPTH"
    
    find "$SC_OUT/amp/MDA/sim*/" -type f -name "*recal.bam" | sort > "$SC_OUT/amp/MDA/downsampled_MDA_bam.txt"
    
    #conda activate SingleCellSim
    bash "$SIGMA_DIR/real_data/subset_sort_index.sh" \
      --input "$SC_OUT/amp/MDA/downsampled_MDA_bam.txt" \
      --cores "$CORES"
    
    echo "scDNAseq simulation complete. Outputs saved to $SC_OUT"
    ;;

  bulk)
    if [ -z "$REF" ]; then
      echo "Error: --ref is required for bulk"
      exit 1
    fi
    echo "=== Running Bulk DNAseq Simulation ==="
    
    # 7. Bulk Simulation
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
