import argparse
import os
import sys
import pickle
import random

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

from BulkSim.bulkFunctions import read_fasta,read_vcf, apply_snps, amplify_mutations, save_fasta


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate bulk sequencing data in FASTA format.")
    parser.add_argument("--reference", type=str, required=True, help="Path to the reference genome in FASTA format.")
    parser.add_argument("--vcf", type=str, required=True, help="Path to the VCF file containing paternal and maternal SNPs.")
    parser.add_argument("--mutations", type=str, required=True, help="Path to the pickle file containing coalescent tree mutations.")
    parser.add_argument("--output", type=str, required=True, help="Output FASTA file for bulk genome.")

    args = parser.parse_args()

    # Read the reference genome
    reference = read_fasta(args.reference)

    # Read the SNPs from the VCF file
    snps = read_vcf(args.vcf)

    # Load the coalescent tree mutations and VAF
    with open(args.mutations, "rb") as f:
        data = pickle.load(f)
    mutations = data["mutations"]
    vaf = data.get("vaf", {})

    # Create paternal and maternal genomes
    paternal_genome = apply_snps(reference, snps, paternal=True)
    maternal_genome = apply_snps(reference, snps, paternal=False)

    # Amplify mutations in paternal and maternal genomes based on VAF
    bulk_genome_paternal = amplify_mutations(paternal_genome, mutations, vaf)
    bulk_genome_maternal = amplify_mutations(maternal_genome, mutations, vaf)

    # Combine paternal and maternal genomes for bulk sequencing
    combined_bulk_genome = {chrom: bulk_genome_paternal[chrom] for chrom in bulk_genome_paternal}
    for chrom in bulk_genome_maternal:
        if chrom in combined_bulk_genome:
            combined_bulk_genome[chrom] = ''.join(
                random.choice([bulk_genome_paternal[chrom][i], bulk_genome_maternal[chrom][i]])
                for i in range(len(bulk_genome_paternal[chrom]))
            )

    # Save the resulting bulk genome as a FASTA file
    save_fasta(combined_bulk_genome, args.output)

    print(f"Bulk genome saved to '{args.output}'.")
