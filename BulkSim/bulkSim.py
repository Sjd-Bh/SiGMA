import argparse
import os
import sys
import pickle

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
bulkSim_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.insert(0, bulkSim_dir)


from BulkSim.bulkFunctions import read_vcf, apply_snps_to_reference, amplify_genomes


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Simulate bulk genome amplification with SNPs and mutations.")
    parser.add_argument("--ref", type=str, required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument("--mat-snp", type=str, required=True, help="Path to the maternal SNPs VCF file.")
    parser.add_argument("--pat-snp", type=str, required=True, help="Path to the paternal SNPs VCF file.")
    parser.add_argument("--coal", type=str, required=True, help="Path to the coalescent data pickle file.")
    parser.add_argument("--output", type=str, required=True, help="Output folder to save the amplified genomes and mutations.")
    args = parser.parse_args()

    # Create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Load the coalescent data
    with open(args.coal, 'rb') as f:
        coalescent_data = pickle.load(f)

    #tree = coalescent_data['tree']
    mutations = coalescent_data['mutations']
    vaf_info = coalescent_data['vaf']

    # Read SNPs from VCF files
    pat_snp_positions = read_vcf(args.pat_snp)
    mat_snp_positions = read_vcf(args.mat_snp)
    
    # Apply SNPs to the reference genome to generate maternal and paternal sequences
    paternal_genome_path = os.path.join(args.output, "paternal_genome.fasta")
    maternal_genome_path = os.path.join(args.output, "maternal_genome.fasta")
    apply_snps_to_reference(args.ref, pat_snp_positions, paternal_genome_path)
    apply_snps_to_reference(args.ref, mat_snp_positions, maternal_genome_path)
    
    # Amplify genomes based on mutations from the coalescent tree and save mutations in a VCF file
    amplified_fasta_path = os.path.join(args.output, "amplified_genomes.fasta")
    amplified_vcf_path = os.path.join(args.output, "amplified_genomes.vcf")
    amplify_genomes(paternal_genome_path, maternal_genome_path, amplified_fasta_path, amplified_vcf_path, mutations, vaf_info)

if __name__ == "__main__":
    main()
