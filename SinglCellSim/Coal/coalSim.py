import pickle
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def coalescent_tree_simulation(n_cells, genome_length, mutation_rate, N=1000, filename="coalescent_tree.pkl"):
    """
    Simulates a coalescent tree for somatic evolution, assigns mutations, 
    propagates them through the tree, saves the tree, and visualizes it.

    Parameters:
    - n_cells: Number of cells to simulate.
    - genome_length: Length of the genome.
    - mutation_rate: Mutation rate per unit branch length.
    - N: Effective population size.
    - filename: File to save the tree simulation (Pickle format).
    """
    # Initialize variables
    active_lineages = list(range(n_cells))
    next_internal_node = n_cells
    branch_lengths = {}
    tree = nx.DiGraph()

    # Simulate tree topology and branch lengths
    while len(active_lineages) > 1:
        # Sample coalescent time
        k = len(active_lineages)
        coalescent_time = np.random.exponential(scale=(4 * N) / (k * (k - 1)))

        # Randomly choose two lineages to merge
        lineage1, lineage2 = np.random.choice(active_lineages, size=2, replace=False)
        active_lineages.remove(lineage1)
        active_lineages.remove(lineage2)

        # Add a new internal node
        tree.add_node(next_internal_node, time=coalescent_time)
        tree.add_edges_from([(next_internal_node, lineage1), (next_internal_node, lineage2)])
        branch_lengths[(next_internal_node, lineage1)] = coalescent_time
        branch_lengths[(next_internal_node, lineage2)] = coalescent_time
        active_lineages.append(next_internal_node)
        next_internal_node += 1

    # Assign mutations and propagate them
    mutations = {}
    node_mutations = {node: set() for node in tree.nodes}
    mutation_positions = set()

    for (parent, child), length in branch_lengths.items():
        num_mutations = np.random.poisson(mutation_rate * length)
        branch_mutations = set()
        for _ in range(num_mutations):
            while True:
                position = np.random.randint(0, genome_length)
                if position not in mutation_positions:
                    mutation_positions.add(position)
                    branch_mutations.add(position)
                    break
        mutations[(parent, child)] = branch_mutations

    # Propagate mutations to all descendants
    def propagate_mutations(node):
        for child in tree.successors(node):
            node_mutations[child].update(node_mutations[node])  # Inherit from parent
            node_mutations[child].update(mutations.get((node, child), set()))  # Add branch-specific mutations
            propagate_mutations(child)

    # Start propagation from the root (last internal node)
    root = next_internal_node - 1
    propagate_mutations(root)

    # Count occurrences of each mutation in leaf nodes only
    # Count the total number of leaf nodes
    leaf_nodes = [node for node in tree.nodes if tree.out_degree(node) == 0]
    
    # Calculate the VAF for each mutation
    vaf_info = {
        pos: sum(pos in node_mutations[leaf] for leaf in leaf_nodes) / len(leaf_nodes)
        for pos in mutation_positions
    }
    # Organize the output data
    output_data = {
        "tree": tree,
        "mutations": mutations,
        "vaf_info": vaf_info,
        "node_mutations": node_mutations  # Store propagated mutations per node
    }

    # Save the output data to a pickle file
    with open(filename, "wb") as f:
        pickle.dump(output_data, f)

    # Visualize and save the tree
    visualize_and_save_tree(tree, filename.replace(".pkl", ".png"))
    print(f"Data saved to {filename}, visualization saved as {filename.replace('.pkl', '.png')}")

def visualize_and_save_tree(tree, image_filename):
    """
    Visualizes and saves the coalescent tree as an image file.

    Parameters:
    - tree: NetworkX DiGraph representation of the tree.
    - image_filename: Path to save the visualization image.
    """
    # Use Graphviz layout for proper tree structure
    pos = nx.drawing.nx_pydot.graphviz_layout(tree, prog="dot")

    # Draw the tree
    plt.figure(figsize=(10, 8))
    nx.draw(
        tree, pos, with_labels=True,
        node_size=1000, node_color="skyblue",
        font_size=10, font_weight="bold",
        arrowsize=15
    )
    plt.title("Coalescent Tree Visualization")
    plt.savefig(image_filename)
    plt.close()

# Example usage
coalescent_tree_simulation(
    n_cells=10,
    genome_length=1000000,
    mutation_rate=1e-3,
    N=1000,
    filename="test_tree.pkl"
)
