#!/usr/bin/env python3
import argparse
import pickle
import numpy as np
import networkx as nx

def coalescent_tree_simulation(n_cells, genome_length, N=1000, filename="coalescent_tree.pkl"):
    """
    Simulates a coalescent tree for somatic evolution, assigns mutations, 
    propagates them through the tree, and saves it as a pickle file.

    Parameters:
    - n_cells (int): Number of cells (leaves) to simulate.
    - genome_length (int): Length of the genome.
    - mutation_rate (float): Mutation rate per unit branch length.
    - N (int): Effective population size.
    - filename (str): Output filename for saving the tree (Pickle format).
    """
    # Initialize variables
    active_lineages = list(range(n_cells))  # leaf nodes are 0, 1, 2, ..., n_cells-1
    next_internal_node = n_cells            # internal nodes start after leaves
    branch_lengths = {}
    tree = nx.DiGraph()

    # Track node times
    node_times = {i: 0.0 for i in active_lineages}

    current_time = 0.0

    # Simulate coalescent process
    while len(active_lineages) > 1:
        k = len(active_lineages)

        # Sample time UNTIL the next coalescent event
        time_interval = np.random.exponential(scale=(4 * N) / (k * (k - 1)))
        current_time += time_interval  # Update absolute time

        # Choose two lineages to merge
        lineage1, lineage2 = np.random.choice(active_lineages, size=2, replace=False)
        active_lineages.remove(lineage1)
        active_lineages.remove(lineage2)

        # Create a new internal node
        tree.add_node(next_internal_node, time=current_time)

        # Add edges
        tree.add_edge(next_internal_node, lineage1)
        tree.add_edge(next_internal_node, lineage2)

        # Update branch lengths (Parent Time - Child Time)
        branch_lengths[(next_internal_node, lineage1)] = current_time - node_times[lineage1]
        branch_lengths[(next_internal_node, lineage2)] = current_time - node_times[lineage2]

        # Update times for new internal node
        node_times[next_internal_node] = current_time

        # Add internal node back to active list
        active_lineages.append(next_internal_node)
        next_internal_node += 1

    # Final root node
    root = active_lineages[0]

    # Save tree to pickle file
    tree_data = {
        "tree": tree,
        "root": root,
        "branch_lengths": branch_lengths,
        "node_times": node_times,
        "n_cells": n_cells,
        "genome_length": genome_length,
        "N": N
    }

    with open(filename, "wb") as f:
        pickle.dump(tree_data, f)

    print(f"[INFO] Coalescent tree successfully saved to '{filename}'")
    print(f"[INFO] Root node: {root}")
    print(f"[INFO] Total nodes: {len(tree.nodes())}")
    print(f"[INFO] Total edges: {len(tree.edges())}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate a coalescent tree and save it as a pickle file.")
    parser.add_argument("--n_cells", type=int, required=True, help="Number of leaf cells to simulate")
    parser.add_argument("--genome_length", type=int, required=True, help="Length of the genome")
    parser.add_argument("--N", type=int, default=1000, help="Effective population size (default: 1000)")
    parser.add_argument("--output", type=str, default="coalescent_tree.pkl", help="Output pickle filename")

    args = parser.parse_args()

    print("simulating tree ...")
    coalescent_tree_simulation(
        n_cells=args.n_cells,
        genome_length=args.genome_length,
        N=args.N,
        filename=args.output
    )
