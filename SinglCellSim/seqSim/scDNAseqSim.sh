#!/bin/bash

# First argument is the reference file
reference=$1

# Shift arguments to process all the provided FASTA files
shift

# Loop through all the FASTA files passed to the script
for fasta in "$@"; do
  # Run each iteration in the background
  (
    # Extract the base name (without directory and extension)
    base_name=$(basename "$fasta" .fasta)

    # Set output folder based on the FASTA file directory
    out_folder=$(dirname "$fasta")

    # Run ART to simulate reads
    art_illumina -l 150 -f 3 -m 200 -s 10 -i "$fasta" -o "$out_folder/$base_name"

    # Run HISAT2 to align reads using the specified reference
    hisat2 -x "$reference" -1 "$out_folder/${base_name}1.fq" -2 "$out_folder/${base_name}2.fq" -S "$out_folder/$base_name.sam"

    # Convert SAM to BAM, sort, and index BAM
    samtools view -bS -o "$out_folder/$base_name.bam" "$out_folder/$base_name.sam"
    samtools sort "$out_folder/$base_name.bam" -o "$out_folder/${base_name}_sort.bam"
    samtools index "$out_folder/${base_name}_sort.bam"

    # Add or replace read groups using Picard
    java -jar ../../picard/picard.jar AddOrReplaceReadGroups \
      I="$out_folder/${base_name}_sort.bam" \
      O="$out_folder/${base_name}_sort_rg.bam" \
      RGID=1 \
      RGLB=library_name \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=simBulk

    # Index the sorted BAM with read groups
    samtools index "$out_folder/${base_name}_sort_rg.bam"

    # Run GATK HaplotypeCaller to call variants using the reference
    gatk HaplotypeCaller -R "$reference" -I "$out_folder/${base_name}_sort_rg.bam" -O "$out_folder/${base_name}_sort_rg.vcf"
  ) &
done

# Wait for all background processes to finish
wait
