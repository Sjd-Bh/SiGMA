from Bio import SeqIO
import argparse
import pickle
from collections import defaultdict
import os
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
bulkSim_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.insert(0, bulkSim_dir)

def verify_mutation_vaf(amplified_fasta, mutations, vaf_info):
    mutation_counts = defaultdict(int)
    total_sequences = 0

    # Read the amplified sequences
    for record in SeqIO.parse(amplified_fasta, "fasta"):
        total_sequences += 1
        sequence = record.seq

        # Count the presence of mutations in this sequence
        for mut in mutations:
            if mut < len(sequence) and sequence[mut] != '-':  # Check that the mutation position is valid
                mutation_counts[mut] += 1

    # Calculate observed VAFs and compare them with expected VAFs
    for mut, expected_vaf in vaf_info.items():
        observed_vaf = mutation_counts[mut-1] / total_sequences
        print(f"Mutation {mut}: Observed VAF = {observed_vaf:.2f}, Expected VAF = {expected_vaf:.2f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify mutation VAFs in amplified genomes.")
    parser.add_argument("--amplified-fasta", type=str, required=True, help="Path to the amplified genomes FASTA file.")
    parser.add_argument("--coal", type=str, required=True, help="Path to the coalescent data pickle file.")
    args = parser.parse_args()

    # Load the coalescent data
    with open(args.coal, 'rb') as f:
        coalescent_data = pickle.load(f)

    mutations = coalescent_data['mutations']
    vaf_info = coalescent_data['vaf_info']

    # Verify mutation VAFs
    verify_mutation_vaf(args.amplified_fasta, mutations, vaf_info)
