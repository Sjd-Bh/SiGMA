import argparse
import glob
import os
import subprocess
import sys

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
bulkSim_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.insert(0, bulkSim_dir)

def run_sccaller(ref, sc_files, vcf_files, output_dir, bulk, core, sccaller_path_py):
    for sc_bam, sc_vcf in zip(sc_files, vcf_files):
        sample_name = os.path.basename(os.path.dirname(sc_bam))
        sample_output_dir = os.path.join(output_dir, sample_name)
        
        # Ensure output subdirectory for each sample exists
        if not os.path.exists(sample_output_dir):
            os.makedirs(sample_output_dir)

        # Run sccaller using conda run in the SingleCellSim environment
        output_vcf = os.path.join(sample_output_dir, "sccaller_{}.vcf".format(sample_name))
        sccaller_command = [
            "conda", "run", "-n", "Sccaller", "python2.7", sccaller_path_py, 
            "-b", sc_bam,
            "-f", ref,
            "-o", output_vcf,
            "--bulk", bulk,
            "-t", "hsnp",
            "-s", sc_vcf,
            "-n", str(core)
        ]
        print("Running sccaller for {}".format(sc_bam))
        subprocess.run(sccaller_command, check=True)

def main():
    parser = argparse.ArgumentParser(description="Run ProSolo and convert VCF to BCF.")
    parser.add_argument("--ref", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--scFiles", required=True, help="Glob pattern for single-cell BAM files.")
    parser.add_argument("--vcf", required=True, help="Glob pattern for VCF files.")
    parser.add_argument("--bulk", required=True, help="Path to bulk BAM file.")
    parser.add_argument("--out", required=True, help="Output directory for ProSolo and BCF files.")
    parser.add_argument("--core", required=int, help="Number of threads")
    parser.add_argument("--sccaller-path-py", required=True, help="Path to the sccaller Python script.")
    args = parser.parse_args()

    # Expand file paths using glob
    sc_files = sorted(glob.glob(args.scFiles))
    vcf_files = sorted(glob.glob(args.vcf))
    
    if len(sc_files) != len(vcf_files):
        print("Error: Number of BAM files does not match the number of VCF files.")
        return

    if not os.path.exists(args.out):
        os.makedirs(args.out)
        
    run_sccaller(args.ref, sc_files, vcf_files, args.out, args.bulk, args.core, args.sccaller_path_py)

if __name__ == "__main__":
    main()
