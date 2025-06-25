import argparse
import subprocess
import os

def extract_homozygous_snps(input_vcf, output_vcf):
    """
    Extract homozygous SNPs and index the output VCF.
    """
    # Step 1: Extract homozygous SNPs
    cmd_extract = [
        "bcftools", "view",
        "-g", "hom",       # Homozygous only
        "-v", "snps",      # SNPs only
        "-o", output_vcf,
        "-O", "v",         # VCF format
        input_vcf
    ]
    print(f"Running: {' '.join(cmd_extract)}")
    subprocess.run(cmd_extract, check=True)

    # Step 2: Compress with bgzip
    print("Compressing VCF with bgzip...")
    subprocess.run(["bgzip", "-f", output_vcf], check=True)

    # Step 3: Index with bcftools
    print("Indexing compressed VCF...")
    subprocess.run(["bcftools", "index", f"{output_vcf}.gz"], check=True)

    print(f"✅ Homozygous SNPs saved and indexed: {output_vcf}.gz")

def main():
    parser = argparse.ArgumentParser(description="Extract homozygous SNPs from bulk VCF using bcftools.")
    parser.add_argument("--input-vcf", required=True, help="Path to input bulk VCF.")
    parser.add_argument("--output-dir", required=True, help="Directory to save output VCF.")
    parser.add_argument("--output-name", default="bulk_homo.vcf", help="Output VCF file name (default: bulk_homo.vcf)")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    output_vcf = os.path.join(args.output_dir, args.output_name)

    extract_homozygous_snps(args.input_vcf, output_vcf)

if __name__ == "__main__":
    main()
