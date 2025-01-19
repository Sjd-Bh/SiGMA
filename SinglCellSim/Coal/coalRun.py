import argparse
import os
import sys
import pickle
from collections import defaultdict  # Can be removed if not used
import networkx as nx

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

from SinglCellSim.Coal.coalSim import coalescent_tree_simulation

def parse_genome_length(length_str):
    """Parse genome length from a string like '100kb' or '1mb'."""
    try:
        if length_str.endswith("kb"):
            return int(float(length_str[:-2]) * 1000)
        elif length_str.endswith("mb"):
            return int(float(length_str[:-2]) * 1e6)
        else:
            return int(length_str)
    except ValueError:
        raise ValueError("Invalid genome length format. Use formats like '100kb', '1mb', or an integer.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a coalescent tree with mutations.")
    parser.add_argument("--SNV-rate", type=float, required=True, help="Mutation rate per base per generation.")
    parser.add_argument("--num-cells", type=int, required=True, help="Number of cells.")
    parser.add_argument("--genome-len", type=str, required=True, help="Length of the genome (e.g., '100kb', '1mb').")
    parser.add_argument("--eff-pop-size", "-N", type=int, default=100, help="Effective population size for the coalescence.")
    parser.add_argument("--output", type=str, required=True, help="Output file to save the tree and mutations.")

    args = parser.parse_args()

    try:
        # Parse genome length
        genome_length = parse_genome_length(args.genome_len)

        # Generate output filenames
        output_tree_filename = args.output if args.output.endswith(".pkl") else f"{args.output}.pkl"
        output_visualization_filename = output_tree_filename.replace(".pkl", ".png")

        # Generate the coalescent tree
        coalescent_tree_simulation(
            n_cells=args.num_cells,
            N=args.eff_pop_size,
            genome_length=genome_length,
            mutation_rate=args.SNV_rate,
            filename=output_tree_filename
        )

        print(f"Coalescence tree and associated data saved to '{output_tree_filename}'.")
        print(f"Tree visualization saved to '{output_visualization_filename}'.")
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
