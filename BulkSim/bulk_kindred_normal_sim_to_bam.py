#!/usr/bin/env python3

import argparse
import os
import subprocess

def run_cmd(cmd):
    """Utility to run a shell command and check for errors."""
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def align_bwa(ref, r1, r2, sam_out, threads):
    """4. Alignment using BWA MEM"""
    cmd = ["bwa", "mem", "-t", str(threads), ref, r1, r2]
    print(f"Running: {' '.join(cmd)} > {sam_out}")
    with open(sam_out, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

def sam_to_sorted_bam(sam_file, sorted_bam, threads):
    """5 & 6. Convert SAM to BAM, Sort and Index"""
    bam_file = sam_file.replace(".sam", ".bam")
    run_cmd(["samtools", "view", "-bS", sam_file, "-o", bam_file])
    run_cmd(["samtools", "sort", "-@", str(threads), bam_file, "-o", sorted_bam])
    run_cmd(["samtools", "index", sorted_bam])
    os.remove(sam_file)
    os.remove(bam_file)

def add_read_groups(bam_in, bam_out, sample_name):
    """7. Add sample name (@RG group)"""
    cmd = [
	#"conda", "run", "-n", "picard",
        "picard", "AddOrReplaceReadGroups",
        f"I={bam_in}",
        f"O={bam_out}",
        "RGID=1",
        "RGLB=lib1",
        "RGPL=illumina",
        "RGPU=unit1",
        f"RGSM={sample_name}"
    ]
    run_cmd(cmd)
    run_cmd(["samtools", "index", bam_out])

def mark_duplicates(bam_in, bam_out, metrics_file):
    """8. Mark duplicates"""
    cmd = [
	#"conda", "run", "-n", "picard",
        "picard", "MarkDuplicates",
        f"I={bam_in}",
        f"O={bam_out}",
        f"M={metrics_file}"
    ]
    run_cmd(cmd)
    run_cmd(["samtools", "index", bam_out])

def bqsr(bam_in, final_bam, ref, known_sites):
    """9 & 10. Base Quality Score Recalibration (BQSR)"""
    recal_table = bam_in.replace(".bam", "_recal_data.table")
    
    # BaseRecalibrator
    cmd1 = [
	#"conda", "run", "-n", "picard",
        "gatk", "BaseRecalibrator",
        "-I", bam_in,
        "-R", ref,
        "--known-sites", known_sites,
        "-O", recal_table
    ]
    run_cmd(cmd1)
    
    # ApplyBQSR (Automatically outputs a sorted/indexed BAM if input is sorted)
    cmd2 = [
	#"conda", "run", "-n", "picard",
        "gatk", "ApplyBQSR",
        "-I", bam_in,
        "-R", ref,
        "--bqsr-recal-file", recal_table,
        "-O", final_bam
    ]
    run_cmd(cmd2)
    
    # Ensure final indexing (Step 10)
    run_cmd(["samtools", "index", final_bam])

def process_sample(prefix, r1, r2, ref, known_sites, output_dir, threads):
    print(f"\n--- Processing Sample: {prefix} ---")
    
    sam = os.path.join(output_dir, f"{prefix}.sam")
    sorted_bam = os.path.join(output_dir, f"{prefix}_sorted.bam")
    rg_bam = os.path.join(output_dir, f"{prefix}_rg.bam")
    md_bam = os.path.join(output_dir, f"{prefix}_md.bam")
    metrics = os.path.join(output_dir, f"{prefix}_md_metrics.txt")
    final_bam = os.path.join(output_dir, f"{prefix}_final_bqsr.bam")

    align_bwa(ref, r1, r2, sam, threads)
    sam_to_sorted_bam(sam, sorted_bam, threads)
    add_read_groups(sorted_bam, rg_bam, prefix)
    mark_duplicates(rg_bam, md_bam, metrics)
    bqsr(md_bam, final_bam, ref, known_sites)
    
    # Cleanup intermediate BAMs to save space (optional)
    for tmp_file in [sorted_bam, sorted_bam + ".bai", rg_bam, rg_bam + ".bai", md_bam, md_bam + ".bai"]:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
            
    print(f"--- Finished {prefix}! Final output: {final_bam} ---")

def main():
    parser = argparse.ArgumentParser(description="Process bulk FASTQ to final BQSR BAM")
    parser.add_argument("--ref", required=True, help="Reference FASTA (must be bwa indexed)")
    parser.add_argument("--known_sites", required=True, help="VCF file for BQSR known sites")
    parser.add_argument("--r1", required=True, help="Read 1 FASTQ")
    parser.add_argument("--r2", required=True, help="Read 2 FASTQ")
    parser.add_argument("--sample_name", required=True, help="Sample name for Read Groups")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    process_sample(
        args.sample_name,
        args.r1,
        args.r2,
        args.ref,
        args.known_sites,
        args.output_dir,
        args.threads
    )

if __name__ == "__main__":
    main()
