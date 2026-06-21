import subprocess
from pathlib import Path
import random
import numpy as np
import argparse
from Bio import SeqIO

def simulate_snps(genome_length, snp_rate=1/1000):
    """Simulate SNP positions in a genome (returning 0-based positions)."""
    expected_snps = genome_length * snp_rate
    total_snps = np.random.poisson(expected_snps)
    snp_positions = np.random.choice(range(genome_length), size=total_snps, replace=False)
    return set(snp_positions)

def get_random_alt_base(ref_base):
    """Return a random base different from the reference."""
    bases = {'A', 'T', 'C', 'G'}
    return random.choice(list(bases - {ref_base.upper()}))

def process_vcf_with_bcftools(vcf_path, ref_path):
    """Sort, norm, bgzip, and index the VCF using bcftools."""
    vcf_path = str(vcf_path)
    sorted_vcf = vcf_path.replace(".vcf", ".sorted.vcf.gz")
    dedup_vcf = vcf_path.replace(".vcf", ".dedup.vcf.gz")

    print(f"Processing {vcf_path} with bcftools...")
    # Sort and compress
    subprocess.run(f"conda run -n picard bcftools sort {vcf_path} -O z -o {sorted_vcf}", shell=True, check=True)
    subprocess.run(f"conda run -n picard bcftools index -t {sorted_vcf}", shell=True, check=True)
    
    # Normalize, left-align, and remove duplicates
    subprocess.run(f"conda run -n picard bcftools norm -f {ref_path} -d both {sorted_vcf} -O z -o {dedup_vcf}", shell=True, check=True)
    subprocess.run(f"conda run -n picard bcftools index -t {dedup_vcf}", shell=True, check=True)
    print(f"Finished: {dedup_vcf}")

def generate_and_save_snps(output_dir, snp_rate, ref_path):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    merged_vcf_path = output_dir / "merged_snps.vcf"
    maternal_vcf_path = output_dir / "maternal_snps.vcf"
    paternal_vcf_path = output_dir / "paternal_snps.vcf"

    # 1. Parse reference to get contigs for the VCF header
    records = list(SeqIO.parse(ref_path, "fasta"))
    contig_headers = "".join([f"##contig=<ID={r.id},length={len(r.seq)}>\n" for r in records])

    base_header = (
        "##fileformat=VCFv4.2\n"
        f"{contig_headers}"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
    )

    with open(merged_vcf_path, "w") as merged_vcf, \
         open(maternal_vcf_path, "w") as mat_vcf, \
         open(paternal_vcf_path, "w") as pat_vcf:

        merged_vcf.write(base_header)
        mat_vcf.write(base_header)
        pat_vcf.write(base_header)

        for record in records:
            chrom = record.id
            chrom_seq = str(record.seq).upper()
            chrom_length = len(chrom_seq)

            maternal_positions = simulate_snps(chrom_length, snp_rate)
            paternal_positions = simulate_snps(chrom_length, snp_rate)

            # Prevent duplicate overlapping variants (fixes GATK multiple variant context error)
            overlaps = maternal_positions.intersection(paternal_positions)
            maternal_positions -= overlaps
            paternal_positions -= overlaps

            all_positions = sorted(list(maternal_positions.union(paternal_positions)))

            for pos in all_positions:
                ref_base = chrom_seq[pos]
                if ref_base not in ['A', 'T', 'C', 'G']:
                    continue # Skip Ns

                alt_base = get_random_alt_base(ref_base)
                
                if pos in maternal_positions:
                    gt = "0|1"
                    vcf_line = f"{chrom}\t{pos+1}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\tGT\t{gt}\n"
                    mat_vcf.write(vcf_line)
                    merged_vcf.write(vcf_line)
                else:
                    gt = "1|0"
                    vcf_line = f"{chrom}\t{pos+1}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\tGT\t{gt}\n"
                    pat_vcf.write(vcf_line)
                    merged_vcf.write(vcf_line)

    print(f"Raw VCFs generated in {output_dir}. Starting bcftools processing...")
    
    # Process the generated VCFs
    process_vcf_with_bcftools(merged_vcf_path, ref_path)
    process_vcf_with_bcftools(maternal_vcf_path, ref_path)
    process_vcf_with_bcftools(paternal_vcf_path, ref_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate SNPs in a genome.")
    parser.add_argument("-o", "--output", required=True, help="Output directory for the SNP files")
    parser.add_argument("-R", "--reference", required=True, help="Path to the reference genome in FASTA format")
    parser.add_argument("--rate", type=float, default=1/1000, help="SNP rate (default: 1 SNP per 1000 bp)")
    args = parser.parse_args()

    generate_and_save_snps(args.output, args.rate, args.reference)
