from pathlib import Path
import random
import pysam
import numpy as np
import argparse

def simulate_snps(genome_length, snp_rate=1/1000):
    """Simulate SNP positions in a genome."""
    expected_snps = genome_length * snp_rate
    total_snps = np.random.poisson(expected_snps)
    snp_positions = np.random.choice(range(1, genome_length + 1), size=total_snps, replace=False)
    return sorted(snp_positions)

def generate_and_save_snps(output_dir, snp_rate, ref_path):
    ref = pysam.FastaFile(ref_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    merged_vcf_path = output_dir / "merged_snps.vcf"
    with open(merged_vcf_path, "w") as merged_vcf:
        merged_vcf.write("##fileformat=VCFv4.2\n")
        merged_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom in ref.references:
            chrom_length = ref.get_reference_length(chrom)
            maternal_positions = simulate_snps(chrom_length, snp_rate)
            paternal_positions = simulate_snps(chrom_length, snp_rate)
            maternal_vcf_path = output_dir / f"{chrom}_maternal_snps.vcf"
            paternal_vcf_path = output_dir / f"{chrom}_paternal_snps.vcf"
            # Save maternal SNPs
            with open(maternal_vcf_path, "w") as maternal_vcf:
                maternal_vcf.write("##fileformat=VCFv4.2\n")
                maternal_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for pos in maternal_positions:
                    ref_base = ref.fetch(chrom, pos, pos + 1)
                    alt_base = random.choice([b for b in "ACGT" if b != ref_base])
                    maternal_vcf.write(f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
                    merged_vcf.write(f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
            # Save paternal SNPs
            with open(paternal_vcf_path, "w") as paternal_vcf:
                paternal_vcf.write("##fileformat=VCFv4.2\n")
                paternal_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for pos in paternal_positions:
                    ref_base = ref.fetch(chrom, pos, pos + 1)
                    alt_base = random.choice([b for b in "ACGT" if b != ref_base])
                    paternal_vcf.write(f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
                    merged_vcf.write(f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
    print(f"Maternal and paternal SNPs saved to {output_dir}, merged VCF: {merged_vcf_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate SNPs in a genome.")
    parser.add_argument("-o", "--output", required=True, help="Output directory for the SNP files")
    parser.add_argument("-R", "--reference", required=True, help="Path to the reference genome in FASTA format")
    parser.add_argument("--rate", type=float, default=1/1000, help="SNP rate (default: 1 SNP per 1000 bp)")
    args = parser.parse_args()
    
    generate_and_save_snps(args.output, args.rate, args.reference)
