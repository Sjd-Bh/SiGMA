#!/bin/bash

# Default values
INPUT_LIST=""
CHR=""
CORES=4 # Default number of parallel jobs
OUTPUT_DIR=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT_LIST="$2"; shift ;;
        --chr) CHR="$2"; shift ;;
        --cores) CORES="$2"; shift ;;
        --output-folder) OUTPUT_DIR="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check if arguments are provided
if [ -z "$INPUT_LIST" ]; then
    echo "Usage: $0 --input <path_to_list.txt> [--chr <chromosome>] [--cores <number>] [--output-folder <path>]"
    echo "Example (subset): $0 --input PTA_list.txt --chr 10 --cores 8 --output-folder /path/to/output"
    echo "Example (whole file): $0 --input PTA_list.txt --cores 8"
    exit 1
fi

if [ ! -f "$INPUT_LIST" ]; then
    echo "Error: Input list file '$INPUT_LIST' not found!"
    exit 1
fi

# Ensure chromosome format is 'chrN' ONLY if a chromosome was provided
if [ -n "$CHR" ] && [[ ! "$CHR" == chr* ]]; then
    CHR="chr${CHR}"
fi

# Determine the output directory
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=$(dirname "$INPUT_LIST")
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Export variables so they can be accessed by the subshells running in parallel
export CHR
export OUTPUT_DIR

if [ -n "$CHR" ]; then
    echo "Starting parallel processing (Subsetting, Sorting, Indexing) for chromosome: $CHR"
else
    echo "Starting parallel processing (Sorting and Indexing whole BAMs)"
fi
echo "Output directory: $OUTPUT_DIR"
echo "Running on $CORES cores..."
echo "----------------------------------------"

# Define the processing function
process_bam() {
    local BAM_PATH="$1"
    local CHR="$2"
    local OUTDIR="$3"

    # Extract just the filename (e.g., SRR8438253.recal.bam)
    local BASE_NAME=$(basename "$BAM_PATH")
    
    # Remove the .bam extension for naming output files
    local PREFIX="${BASE_NAME%.bam}"

    # Check if the input file exists using the EXACT path from the text file
    if [ ! -f "$BAM_PATH" ]; then
        echo "Warning: Input BAM '$BAM_PATH' not found. Skipping..."
        return
    fi

    echo "Processing $BASE_NAME..."
    
    if [ -n "$CHR" ]; then
        # CHR is provided: Subset, Sort, Index
        local OUT_BAM="${OUTDIR}/${PREFIX}_${CHR}.bam"
        local SORTED_BAM="${OUTDIR}/${PREFIX}_${CHR}.sorted.bam"

        samtools view -b "$BAM_PATH" "$CHR" > "$OUT_BAM"
        samtools sort "$OUT_BAM" -o "$SORTED_BAM"
        samtools index "$SORTED_BAM"
        
        # Clean up unsorted subset
        rm "$OUT_BAM"
    else
        # No CHR provided: Just Sort and Index
        local SORTED_BAM="${OUTDIR}/${PREFIX}.sorted.bam"

        samtools sort "$BAM_PATH" -o "$SORTED_BAM"
        samtools index "$SORTED_BAM"
    fi
}
export -f process_bam

# Read the list, ignore empty lines, and pipe to xargs for parallel execution
grep -v '^$' "$INPUT_LIST" | xargs -P "$CORES" -I {} bash -c 'process_bam "{}" "$CHR" "$OUTPUT_DIR"'

echo "----------------------------------------"
echo "All samples processed successfully."
