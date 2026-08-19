#!/bin/bash

# ==============================================================================
# Script 2: Picard Downsample and GATK BQSR Pipeline
# Sequential execution per sample to prevent OutOfMemory errors.
# ==============================================================================

# Initialize variables
INPUT_LIST=""
PICARD=""
REF=""
MILLS=""
DBSNP=""
CSV_FILE=""
TARGET_DEPTH=""
MAX_JOBS=4  # CHANGED from 5 to 1 to prevent RAM exhaustion

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT_LIST="$2"; shift ;;
        --picard) PICARD="$2"; shift ;;
        --ref) REF="$2"; shift ;;
        --mills) MILLS="$2"; shift ;;
        --dbsnp) DBSNP="$2"; shift ;;
        --csv) CSV_FILE="$2"; shift ;;
        --target-depth) TARGET_DEPTH="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate inputs (MILLS is now excluded from the mandatory check)
if [ -z "$INPUT_LIST" ] || [ -z "$PICARD" ] || [ -z "$REF" ] || [ -z "$DBSNP" ] || [ -z "$CSV_FILE" ] || [ -z "$TARGET_DEPTH" ]; then
    echo "Error: Missing arguments."
    echo "Usage: ./downsamp_BQ.sh --input <bam.txt> --picard <picard.jar> --ref <ref.fa> --dbsnp <dbsnp.vcf> [--mills <mills.vcf>] --csv <combined_depths.csv> --target-depth <min_depth>"
    exit 1
fi

if [ ! -f "$INPUT_LIST" ]; then echo "Error: Input file '$INPUT_LIST' not found!"; exit 1; fi
if [ ! -f "$CSV_FILE" ]; then echo "Error: CSV file '$CSV_FILE' not found!"; exit 1; fi

# Define the processing function
process_pipeline() {
    local input_bam="$1"
    local picard_jar="$2"
    local ref_fa="$3"
    local mills_vcf="$4"
    local dbsnp_vcf="$5"
    local csv_file="$6"
    local target_depth="$7"

    local input_dir=$(dirname "$input_bam")
    local filename=$(basename "$input_bam")
    local sample_name="${filename%%.*}"

    local downsampled_bam="${input_dir}/${sample_name}.downsampled.bam"
    local recal_table="${input_dir}/${sample_name}.recal.table"
    local recal_bam="${input_dir}/${sample_name}.recal.bam"

    echo "--------------------------------------------------------"
    echo "Starting pipeline for sample: $sample_name"

    # 1. Get sample depth from CSV and calculate probability P
    # Assumes CSV format: Sample_Name,Mean_Depth,BAM_Path
    local sample_depth=$(awk -F, -v s="$sample_name" '$1 == s {print $2}' "$csv_file" | head -n 1)

    if [ -z "$sample_depth" ]; then
        echo "ERROR: Depth for $sample_name not found in $csv_file. Skipping."
        return 1
    fi

    # Calculate P = Target / Sample_Depth. If P > 1, set to 1.
    local prob=$(awk -v t="$target_depth" -v d="$sample_depth" 'BEGIN { p = t/d; if(p>1) p=1; print p }')

    echo "Sample Depth: $sample_depth | Target Depth: $target_depth | Downsample Prob: $prob"

    # 2. Downsample BAM
    echo "Running Picard DownsampleSam for $sample_name..."
    java -Xmx16G -jar "$picard_jar" DownsampleSam \
        I="$input_bam" \
        O="$downsampled_bam" \
        P="$prob" \
        STRATEGY=ConstantMemory \
        MAX_RECORDS_IN_RAM=150000 \
        VALIDATION_STRINGENCY=SILENT

    if [ $? -ne 0 ]; then echo "ERROR: DownsampleSam failed for $sample_name"; return 1; fi

    # 3. GATK BaseRecalibrator
    echo "Running GATK BaseRecalibrator for $sample_name..."
    
    # Conditionally build the known-sites array
    local KNOWN_SITES=("--known-sites" "$dbsnp_vcf")
    if [ -n "$mills_vcf" ]; then
        KNOWN_SITES+=("--known-sites" "$mills_vcf")
    fi

    gatk --java-options "-Xmx16G" BaseRecalibrator \
        -R "$ref_fa" \
        -I "$downsampled_bam" \
        "${KNOWN_SITES[@]}" \
        -O "$recal_table"

    if [ $? -ne 0 ]; then echo "ERROR: BaseRecalibrator failed for $sample_name"; return 1; fi

    # 4. GATK ApplyBQSR
    echo "Running GATK ApplyBQSR for $sample_name..."
    gatk --java-options "-Xmx16G" ApplyBQSR \
        -R "$ref_fa" \
        -I "$downsampled_bam" \
        --bqsr-recal-file "$recal_table" \
        -O "$recal_bam"

    if [ $? -eq 0 ]; then
        echo "SUCCESS: Completed pipeline for $sample_name."
    else
        echo "ERROR: ApplyBQSR failed for $sample_name."
        return 1
    fi
}

# Export function and variables so xargs can use them
export -f process_pipeline

echo "Starting Downsample and BQSR pipeline with $MAX_JOBS concurrent jobs..."

# Run xargs
cat "$INPUT_LIST" | xargs -n 1 -P "$MAX_JOBS" -I {} bash -c 'process_pipeline "{}" "'"$PICARD"'" "'"$REF"'" "'"$MILLS"'" "'"$DBSNP"'" "'"$CSV_FILE"'" "'"$TARGET_DEPTH"'"'

echo "All jobs completed."
