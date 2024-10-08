import os
import argparse
import subprocess
from multiprocessing import Pool
import sys

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

def process_fasta(reference, fasta):
    base_name = os.path.basename(fasta).replace('.fasta', '')
    out_folder = os.path.dirname(fasta)
    
    # Run ART to simulate reads
    art_cmd = [
        'art_illumina', '-l', '150', '-f', '3', '-m', '200', '-s', '10',
        '-i', fasta, '-o', f"{out_folder}/{base_name}"
    ]
    subprocess.run(art_cmd)
    
    # Run HISAT2 to align reads
    hisat2_cmd = [
        'hisat2', '-x', reference, '-1', f"{out_folder}/{base_name}1.fq",
        '-2', f"{out_folder}/{base_name}2.fq", '-S', f"{out_folder}/{base_name}.sam"
    ]
    subprocess.run(hisat2_cmd)
    
    # Convert SAM to BAM, sort, and index
    samtools_view_cmd = ['samtools', 'view', '-bS', '-o', f"{out_folder}/{base_name}.bam", f"{out_folder}/{base_name}.sam"]
    subprocess.run(samtools_view_cmd)
    
    samtools_sort_cmd = ['samtools', 'sort', f"{out_folder}/{base_name}.bam", '-o', f"{out_folder}/{base_name}_sort.bam"]
    subprocess.run(samtools_sort_cmd)
    
    samtools_index_cmd = ['samtools', 'index', f"{out_folder}/{base_name}_sort.bam"]
    subprocess.run(samtools_index_cmd)
    
    # Add or replace read groups using Picard
    picard_cmd = [
        'java', '-jar', '../../picard/picard.jar', 'AddOrReplaceReadGroups',
        'I=' + f"{out_folder}/{base_name}_sort.bam",
        'O=' + f"{out_folder}/{base_name}_sort_rg.bam",
        'RGID=1', 'RGLB=library_name', 'RGPL=illumina', 'RGPU=unit1', 'RGSM=simBulk'
    ]
    subprocess.run(picard_cmd)
    
    # Index sorted BAM with read groups
    samtools_index_rg_cmd = ['samtools', 'index', f"{out_folder}/{base_name}_sort_rg.bam"]
    subprocess.run(samtools_index_rg_cmd)
    
    # Run GATK HaplotypeCaller
    gatk_cmd = [
        'gatk', 'HaplotypeCaller', '-R', reference, '-I', f"{out_folder}/{base_name}_sort_rg.bam",
        '-O', f"{out_folder}/{base_name}_sort_rg.vcf"
    ]
    subprocess.run(gatk_cmd)

def main():
    parser = argparse.ArgumentParser(description="Run ART, HISAT2, and GATK HaplotypeCaller for multiple FASTA files")
    parser.add_argument('--ref', required=True, help='Path to the reference file (FASTA)')
    parser.add_argument('--scFiles', required=True, nargs='+', help='Paths to the single-cell FASTA files')
    parser.add_argument('--cores', type=int, default=1, help='Number of cores to use (default: 1)')
    
    args = parser.parse_args()
    
    # Use multiprocessing to run the process_fasta function with the specified number of cores
    with Pool(args.cores) as pool:
        pool.starmap(process_fasta, [(args.ref, fasta) for fasta in args.scFiles])

if __name__ == '__main__':
    main()
