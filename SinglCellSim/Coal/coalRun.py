import argparse
import os
import sys
import pickle
from collections import defaultdict

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

# Import functions and classes
from SinglCellSim.Coal.coalescnet_new_instead import TreeNode
from SinglCellSim.Coal.coalescnet_new_instead import CoalescentTree
from SinglCellSim.Coal.coalescnet_new_instead import plot_coalescent_tree

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a coalescent tree with mutations, CNVs, and selective sweeps.")
    parser.add_argument("--SNV-rate", type=float, required=True, help="Mutation rate per base per generation.")
    parser.add_argument("--num-cells", type=int, required=True, help="Number of cells.")
    parser.add_argument("--genome-len", type=str, required=True, help="Length of the genome (e.g., '100kb', '1mb').")
    parser.add_argument("--CNV-rate", type=float, default=0.01, help="Rate of CNV events.")
    parser.add_argument("--CNV-size-mean", type=int, default=1000, help="Mean size of CNVs.")
    parser.add_argument("--CNV-size-std", type=int, default=500, help="Standard deviation of CNV sizes.")
    parser.add_argument("--sweep-strength", type=float, default=0.05, help="Probability of selective sweep occurring.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--eff-pop-size", "-N", type=int, default=100, help="Effective population size for the coalescence.")
    parser.add_argument("--output", type=str, required=True, help="Output file to save the tree and mutations.")
    # parser.add_argument("--plot-output", type=str, required=True, help="Output path to save the tree visualization.")
    
    args = parser.parse_args()

    # Parse genome length
    if args.genome_len.endswith("kb"):
        genome_length = int(float(args.genome_len[:-2]) * 1000)
    elif args.genome_len.endswith("mb"):
        genome_length = int(float(args.genome_len[:-2]) * 1e6)
    else:
        try:
            genome_length = int(args.genome_len)
        except ValueError:
            print("Error: Invalid genome length format. Use formats like '100kb', '1mb', or a numeric value.")
            sys.exit(1)

    # Generate the coalescent tree
    tree_simulator = CoalescentTree(
        num_cells=args.num_cells,
        N=args.eff_pop_size,
        seed=args.seed,
        genome_length=genome_length,
        mutation_rate=args.SNV_rate,
        CNV_rate=args.CNV_rate,
        CNV_size_mean=args.CNV_size_mean,
        CNV_size_std=args.CNV_size_std,
        sweep_strength=args.sweep_strength
    )
    tree_root = tree_simulator.tree

    # Collect all leaf nodes
    leaf_nodes = []
    tree_simulator.collect_leaf_nodes(tree_root, leaf_nodes)

    # Count the occurrences of each mutation in leaf nodes only
    leaf_mutation_counts = defaultdict(int)
    for leaf in leaf_nodes:
        for mutation in leaf.mutations:
            leaf_mutation_counts[mutation] += 1

    # Calculate VAFs based on leaf nodes only
    vaf_info = {mutation: count / len(leaf_nodes) for mutation, count in leaf_mutation_counts.items()}

    # Collect mutations for each node
    node_mutations = {}
    tree_simulator.collect_node_mutations(tree_root, node_mutations)

    # Save the tree, mutations, VAFs, and node mutations to the output file
    output_data = {
        "tree": tree_root,
        "mutations": list(leaf_mutation_counts.keys()),
        "vaf_info": vaf_info,
        "node_mutations": node_mutations,
    }

    try:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, "wb") as f:
            pickle.dump(output_data, f)
        print(f"Coalescent tree, mutations, VAFs, and node mutations saved to '{args.output}'.")
    except Exception as e:
        print(f"Error saving output: {e}")
        sys.exit(1)

    # Plot and save the tree visualization
    # try:
    #     os.makedirs(os.path.dirname(args.plot_output), exist_ok=True)
    #     plot_coalescent_tree(tree_root, args.plot_output)
    #     print(f"Tree visualization saved to '{args.plot_output}'.")
    # except Exception as e:
    #     print(f"Error saving plot: {e}")
    #     sys.exit(1)
