###########################################
import argparse
import os
import sys
import pickle
from collections import defaultdict

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

from SinglCellSim.Coal.caolFunctions import make_coalescence_tree, collect_mutations, propagate_mutations, collect_mutation_counts, collect_node_mutations, collect_leaf_nodes

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a coalescence tree with mutations.")
    parser.add_argument("--SNV-rate", type=float, required=True, help="Mutation rate per base per generation.")
    parser.add_argument("--num-cells", type=int, required=True, help="Number of cells.")
    parser.add_argument("--genome-len", type=str, required=True, help="Length of the genome (e.g., '100kb', '1mb').")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--eff-pop-size", "-N", type=int, default=100, help="Effective population size for the coalescence.")
    parser.add_argument("--output", type=str, required=True, help="Output file to save the tree and mutations.")
    
    args = parser.parse_args()

    # Parse genome length
    if args.genome_len.endswith("kb"):
        genome_length = int(float(args.genome_len[:-2]) * 1000)
    elif args.genome_len.endswith("mb"):
        genome_length = int(float(args.genome_len[:-2]) * 1e6)
    else:
        genome_length = int(args.genome_len)

    # Generate the coalescence tree
    tree_root = make_coalescence_tree(args.num_cells, N=args.eff_pop_size, seed=args.seed, genome_length=genome_length, mutation_rate=args.SNV_rate)

    # Collect all mutations
    mutations = set()
    collect_mutations(tree_root, mutations)

    # Collect all leaf nodes
    leaf_nodes = []
    collect_leaf_nodes(tree_root, leaf_nodes)
    
    # Count the occurrences of each mutation in leaf nodes only
    leaf_mutation_counts = defaultdict(int)
    for leaf in leaf_nodes:
        for mutation in leaf.mutations:
            leaf_mutation_counts[mutation] += 1
    
    # Calculate VAFs based on leaf nodes only
    vaf_info = {mutation: count / len(leaf_nodes) for mutation, count in leaf_mutation_counts.items()}

    # Collect mutations for each node
    node_mutations = {}
    collect_node_mutations(tree_root, node_mutations)

    # Save the tree, mutations, VAFs, and node mutations to the output file
    output_data = {
        "tree": tree_root,
        "mutations": list(mutations),
        "vaf_info": vaf_info,
        "node_mutations": node_mutations
    }

    with open(args.output, "wb") as f:
        pickle.dump(output_data, f)

    print(f"Coalescence tree, mutations, VAFs, and node mutations saved to '{args.output}'.")
