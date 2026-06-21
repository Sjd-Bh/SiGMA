import os
import argparse
import subprocess
from multiprocessing import Pool
import sys

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)


def process_fasta(reference_index, fasta):
    """
    Processes a single FASTA file through the simulation and alignment pipeline.
    1. Simulates Illumina reads using ART.
    2. Aligns reads using BWA-MEM.
    3. Converts SAM to BAM, sorts, and adds read groups.
    4. Indexes the final BAM file.
    """
    try:
        base_name = os.path.basename(fasta).replace('.fasta', '').replace('.fa', '')
        out_folder = os.path.dirname(fasta)
        parent_folder_name = os.path.basename(out_folder)
        
        print(f"--- Processing {base_name} in {parent_folder_name} ---")

        # Define file paths
        read1_fq = f"{out_folder}/{base_name}1.fq"
        read2_fq = f"{out_folder}/{base_name}2.fq"
        output_bam = f"{out_folder}/{base_name}.bam"
        sorted_bam = f"{out_folder}/{base_name}_sort.bam"
        
        # FINAL BAM NAME UPDATED HERE:
        final_bam = f"{out_folder}/{parent_folder_name}_{base_name}.bam"

        # 1. Run ART to simulate reads
        print(f"[1/5] Simulating reads for {base_name} with ART...")
        art_cmd = [
            'art_illumina', '-ss', 'HS25', '-l', '150', '-f', '1', '-m', '200', '-s', '10',
            '-i', fasta, '-o', f"{out_folder}/{base_name}"
        ]
        subprocess.run(art_cmd, check=True, capture_output=True, text=True)

        # 2. Run BWA-MEM and pipe to samtools to create a BAM file
        print(f"[2/5] Aligning reads with BWA and converting to BAM...")
        bwa_cmd = ['bwa', 'mem', reference_index, read1_fq, read2_fq]
        samtools_view_cmd = ['samtools', 'view', '-bS', '-o', output_bam, '-']

        # Start the BWA process
        bwa_process = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
        # Start the samtools process, taking input from BWA's output
        samtools_process = subprocess.Popen(samtools_view_cmd, stdin=bwa_process.stdout)
        
        # Allow bwa_process to receive a SIGPIPE if samtools_process exits
        bwa_process.stdout.close()
        
        # Wait for samtools to finish
        samtools_process.wait()
        
        # Check for errors
        if samtools_process.returncode != 0:
            raise subprocess.CalledProcessError(samtools_process.returncode, samtools_view_cmd)
        if bwa_process.wait() != 0:
            raise subprocess.CalledProcessError(bwa_process.returncode, bwa_cmd)


        # 3. Sort BAM
        print(f"[3/5] Sorting BAM file...")
        samtools_sort_cmd = ['samtools', 'sort', output_bam, '-o', sorted_bam]
        subprocess.run(samtools_sort_cmd, check=True)

        # 4. Add or replace read groups using GATK
        print(f"[4/5] Adding read groups with GATK...")
        gatk_rg_cmd = [
            'gatk', 'AddOrReplaceReadGroups',
            f'I={sorted_bam}',
            f'O={final_bam}',
            'RGID=1', 'RGLB=lib1', 'RGPL=illumina', 'RGPU=unit1', f'RGSM={parent_folder_name}_{base_name}'
        ]
        subprocess.run(gatk_rg_cmd, check=True, capture_output=True, text=True)

        # 5. Index the final BAM file
        print(f"[5/5] Indexing final BAM file...")
        subprocess.run(['samtools', 'index', final_bam], check=True)

        print(f"--- Successfully finished processing {base_name} -> {os.path.basename(final_bam)} ---")

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] A command failed while processing {fasta}:", file=sys.stderr)
        print(f"  Command: {' '.join(e.cmd)}", file=sys.stderr)
        if e.stdout:
            print(f"  Stdout: {e.stdout}", file=sys.stderr)
        if e.stderr:
            print(f"  Stderr: {e.stderr}", file=sys.stderr)

    except Exception as e:
        print(f"[ERROR] An unexpected error occurred while processing {fasta}: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Run ART read simulation, BWA alignment, and GATK post-processing."
    )
    parser.add_argument(
        '--ref', required=True, 
        help='Path to the BWA reference index (the FASTA file you used for `bwa index`).'
    )
    parser.add_argument(
        '--scFiles', required=True, nargs='+', 
        help='One or more single-cell FASTA files to process.'
    )
    parser.add_argument(
        '--cores', type=int, default=1, 
        help='Number of cores for parallel processing (default: 1).'
    )

    args = parser.parse_args()

    # Before starting, check if the BWA index files exist.
    # BWA creates several files, we'll just check for a common one like '.bwt'
    if not os.path.exists(f"{args.ref}.bwt"):
        print(f"Error: BWA index file '{args.ref}.bwt' not found.", file=sys.stderr)
        print("Please ensure you have indexed your reference genome with 'bwa index <ref.fasta>'", file=sys.stderr)
        sys.exit(1)


    print(f"Starting pipeline for {len(args.scFiles)} file(s) using {args.cores} core(s).")
    
    # Create a list of arguments for each parallel task
    tasks = [(args.ref, fasta) for fasta in args.scFiles]

    # Run in parallel
    with Pool(args.cores) as pool:
        pool.starmap(process_fasta, tasks)
        
    print("All processing complete.")


if __name__ == '__main__':
    main()
