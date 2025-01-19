from pathlib import Path
import random
import pysam
import numpy as np

def simulate_snps(genome_length, snp_rate=1/1000):
    """
    Simulate SNP positions in a genome.
    
    Parameters:
    - genome_length (int): Length of the genome in base pairs.
    - snp_rate (float): Expected SNP rate (SNPs per base pair). Default is 1 SNP per 1000 bp.
    
    Returns:
    - snp_positions (list): List of SNP positions in the genome (1-based indexing).
    """
    # Calculate the expected number of SNPs
    expected_snps = genome_length * snp_rate

    # Draw the total number of SNPs from a Poisson distribution
    total_snps = np.random.poisson(expected_snps)

    # Randomly assign SNP positions (1-based indexing)
    snp_positions = np.random.choice(range(1, genome_length + 1), size=total_snps, replace=False)

    return sorted(snp_positions)


def generate_and_save_snps(output_dir, num_snps, ref):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    merged_vcf_path = output_dir / "merged_snps.vcf"
    with open(merged_vcf_path, "w") as merged_vcf:
        # Write VCF header
        merged_vcf.write("##fileformat=VCFv4.2\n")
        merged_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom in ref.references:
            chrom_length = ref.get_reference_length(chrom)
            maternal_positions = simulate_snps(chrom_length)
            paternal_positions = simulate_snps(chrom_length)
            maternal_vcf_path = output_dir / f"{chrom}_maternal_snps.vcf"
            paternal_vcf_path = output_dir / f"{chrom}_paternal_snps.vcf"
            # Save maternal SNPs
            with open(maternal_vcf_path, "w") as maternal_vcf:
                maternal_vcf.write("##fileformat=VCFv4.2\n")
                maternal_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for pos in maternal_positions:
                    ref_base = ref.fetch(chrom, pos, pos + 1)
                    alt_base = random.choice([b for b in "ACGT" if b != ref_base])
                    maternal_vcf.write(f"{chrom}\t{pos + 1}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
                    # Write to merged VCF
                    merged_vcf.write(f"{chrom}\t{pos + 1}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
            # Save paternal SNPs
            with open(paternal_vcf_path, "w") as paternal_vcf:
                paternal_vcf.write("##fileformat=VCFv4.2\n")
                paternal_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for pos in paternal_positions:
                    ref_base = ref.fetch(chrom, pos, pos + 1)
                    alt_base = random.choice([b for b in "ACGT" if b != ref_base])
                    paternal_vcf.write(f"{chrom}\t{pos + 1}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
                    # Write to merged VCF
                    merged_vcf.write(f"{chrom}\t{pos + 1}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
    print(f"Maternal and paternal SNPs saved separately, and merged VCF saved to {merged_vcf_path}")

# Load the reference genome
reference_file = "reference_sequence_600kb.fasta"
reference = pysam.FastaFile(reference_file)

# Call the function with the correct arguments
generate_and_save_snps("output_snps", 600, reference)