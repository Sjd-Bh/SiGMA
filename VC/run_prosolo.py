import argparse
import glob
import os
import subprocess

def run_prosolo(ref, sc_files, vcf_files, output_dir):
    for sc_bam, sc_vcf in zip(sc_files, vcf_files):
        # Get sample name from the path structure (e.g., `200/sampleName`)
        sample_name = os.path.basename(os.path.dirname(os.path.dirname(sc_bam)))
        sample_output_dir = os.path.join(output_dir, sample_name)
        
        # Ensure output subdirectory for each sample exists
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Convert input VCF to BCF
        input_bcf = os.path.join(sample_output_dir, f"{sample_name}_input.bcf")
        print(f"Converting input VCF to BCF for {sc_vcf}")
        bcftools_command = ["bcftools", "view", "-O", "b", "-o", input_bcf, sc_vcf]
        subprocess.run(bcftools_command, check=True)

        # Run ProSolo using the BCF candidate file
        output_bcf = os.path.join(sample_output_dir, f"prosolo_{sample_name}.bcf")
        prosolo_command = [
            "prosolo", "single-cell-bulk", "--omit-indels",
            "--sc-isize-mean", "150", "--sc-isize-sd", "10",
            sc_bam,  # Single-cell BAM
            os.path.join(output_dir, "bulk_genome_sort_rg.bam"),  # Replace with correct bulk BAM path if different
            ref,
            "--candidates", input_bcf,
            "--output", output_bcf
        ]
        print(f"Running ProSolo for {sc_bam}")
        subprocess.run(prosolo_command, check=True)
        
        # Convert ProSolo output BCF to VCF
        output_vcf = os.path.splitext(output_bcf)[0] + ".vcf"
        print(f"Converting ProSolo output BCF to VCF for {output_bcf}")
        bcftools_command = ["bcftools", "view", "-O", "v", "-o", output_vcf, output_bcf]
        subprocess.run(bcftools_command, check=True)

def main():
    parser = argparse.ArgumentParser(description="Run ProSolo and convert VCF to BCF.")
    parser.add_argument("--ref", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--scFiles", required=True, help="Glob pattern for single-cell BAM files.")
    parser.add_argument("--vcf", required=True, help="Glob pattern for VCF files.")
    parser.add_argument("--out", required=True, help="Output directory for ProSolo and BCF files.")
    args = parser.parse_args()

    # Collect BAM and VCF files using glob
    sc_files = sorted(glob.glob(args.scFiles))
    vcf_files = sorted(glob.glob(args.vcf))
    
    if len(sc_files) != len(vcf_files):
        print("Error: Number of BAM files does not match the number of VCF files.")
        return

    os.makedirs(args.out, exist_ok=True)
    run_prosolo(args.ref, sc_files, vcf_files, args.out)

if __name__ == "__main__":
    main()
