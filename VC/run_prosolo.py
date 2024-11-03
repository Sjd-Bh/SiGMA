import argparse
import glob
import os
import subprocess
import sys

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
bulkSim_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.insert(0, bulkSim_dir)

def run_prosolo(ref, sc_files, vcf_files, output_dir, bulk):
    for sc_bam, sc_vcf in zip(sc_files, vcf_files):
        sample_name = os.path.basename(os.path.dirname(sc_bam))
        sample_output_dir = os.path.join(output_dir, sample_name)
        
        # Ensure output subdirectory for each sample exists
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Convert input VCF to BCF using bcftools environment
        input_bcf = os.path.join(sample_output_dir, f"{sample_name}_input.bcf")
        print(f"Converting input VCF to BCF for {sc_vcf}")
        bcftools_command = f"source activate bcftools && bcftools view -O b -o {input_bcf} {sc_vcf}"
        subprocess.run(["bash", "-c", bcftools_command], check=True)

        # Run ProSolo in SingleCellSim environment
        output_bcf = os.path.join(sample_output_dir, f"prosolo_{sample_name}.bcf")
        prosolo_command = (
            f"source activate SingleCellSim && prosolo single-cell-bulk --omit-indels "
            f"{sc_bam} {bulk} {ref} --candidates {input_bcf} --output {output_bcf} "
            "--sc-isize-mean 150 --sc-isize-sd 10"
        )
        print(f"Running ProSolo for {sc_bam}")
        subprocess.run(["bash", "-c", prosolo_command], check=True)
        
        # Convert ProSolo output BCF to VCF using bcftools environment
        output_vcf = os.path.splitext(output_bcf)[0] + ".vcf"
        print(f"Converting ProSolo output BCF to VCF for {output_bcf}")
        bcftools_command = f"source activate bcftools && bcftools view -O v -o {output_vcf} {output_bcf}"
        subprocess.run(["bash", "-c", bcftools_command], check=True)

def main():
    parser = argparse.ArgumentParser(description="Run ProSolo and convert VCF to BCF.")
    parser.add_argument("--ref", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--scFiles", required=True, help="Glob pattern for single-cell BAM files.")
    parser.add_argument("--vcf", required=True, help="Glob pattern for VCF files.")
    parser.add_argument("--bulk", required=True, help="Glob pattern for bulk BAM file")
    parser.add_argument("--out", required=True, help="Output directory for ProSolo and BCF files.")
    args = parser.parse_args()

    # Expand file paths using glob
    sc_files = sorted(glob.glob(args.scFiles))
    vcf_files = sorted(glob.glob(args.vcf))
    
    if len(sc_files) != len(vcf_files):
        print("Error: Number of BAM files does not match the number of VCF files.")
        return

    os.makedirs(args.out, exist_ok=True)
    run_prosolo(args.ref, sc_files, vcf_files, args.out, args.bulk)

if __name__ == "__main__":
    main()
