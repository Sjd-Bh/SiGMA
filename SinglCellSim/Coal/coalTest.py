import random
import math
import matplotlib.pyplot as plt
import numpy as np

class TreeNode:
    def __init__(self, name):
        self.left = None
        self.right = None
        self.time = 0.0
        self.name = name
        self.mutations = []

def display_tree(node, level=0):
    if node is not None:
        print(' ' * (level + 4) + node.name + f" (time: {node.time}, mutations: {node.mutations})")
        if node.left or node.right:
            display_tree(node.left, level + 1)
            display_tree(node.right, level + 1)

def propagate_mutations(node, inherited_mutations):
    if node is None:
        return
    # Add inherited mutations without modifying the original list
    node.mutations = inherited_mutations + node.mutations
    # Propagate mutations to children
    propagate_mutations(node.left, node.mutations)
    propagate_mutations(node.right, node.mutations)

def make_coalescence_tree(num_cells, N, seed, genome_length, mutation_rate):
    random.seed(seed)
    np.random.seed(seed)

    active_nodes = [TreeNode(f"cell_{i}") for i in range(num_cells)]
    for node in active_nodes:
        node.time = 0

    next_available_index = num_cells

    while len(active_nodes) > 1:
        total_rate = len(active_nodes) * (len(active_nodes) - 1) / (2 * N)
        time_to_next_event = -math.log(random.random()) / total_rate

        # Accumulate mutations based on the time interval and mutation rate
        for node in active_nodes:
            branch_length = time_to_next_event
            num_mutations = np.random.poisson(branch_length * mutation_rate * genome_length)
            new_mutations = [random.randint(0, genome_length - 1) for _ in range(num_mutations)]
            # Use a new list for mutations to avoid modifying the original mutations list
            node.mutations = node.mutations + new_mutations

        # Choose two nodes to coalesce
        chosen_pair = random.sample(active_nodes, 2)

        # Create a new ancestor node
        new_node = TreeNode(f"anc_{next_available_index}")
        new_node.left, new_node.right = chosen_pair
        new_node.time = chosen_pair[0].time + time_to_next_event

        # No need to inherit mutations here, it will be handled by propagate_mutations

        # Update active nodes
        active_nodes.remove(chosen_pair[0])
        active_nodes.remove(chosen_pair[1])
        active_nodes.append(new_node)

        next_available_index += 1

    # Propagate mutations from the root to all descendant leaves
    propagate_mutations(active_nodes[0], [])

    return active_nodes[0]

def plot_tree(node, x, y, dx):
    if node is not None:
        plt.text(x, y, node.name, ha='center', va='center', bbox=dict(facecolor='white', edgecolor='black'))
        if node.left:
            plt.plot([x, x - dx], [y - 1, y - 5], 'k-')
            plot_tree(node.left, x - dx, y - 6, dx / 2)
        if node.right:
            plt.plot([x, x + dx], [y - 1, y - 5], 'k-')
            plot_tree(node.right, x + dx, y - 6, dx / 2)

# Parameters
num_cells = 10
N = 100
seed = 42
genome_length = 1000000  # Length of the genome
mutation_rate = 1e-8     # Mutation rate per base per generation

# Generate tree and display it
tree_root = make_coalescence_tree(num_cells, N, seed, genome_length, mutation_rate)
display_tree(tree_root)

# Plot the tree
plt.figure(figsize=(10, 10))
plot_tree(tree_root, 0, 0, 20)
plt.gca().axis('off')
plt.show()

# Collect mutations from all cells
def collect_mutations(node, cell_mutations):
    """Collect mutations from all leaf nodes (cells) in the tree."""
    if node.left is None and node.right is None:
        # Leaf node, representing a single cell
        cell_mutations[node.name] = node.mutations
    else:
        # Recursively collect mutations from child nodes
        if node.left:
            collect_mutations(node.left, cell_mutations)
        if node.right:
            collect_mutations(node.right, cell_mutations)

def calculate_vaf(cell_mutations, num_cells):
    """Calculate VAF for each mutation across all cells."""
    mutation_counts = {}
    
    # Count the number of cells containing each mutation
    for mutations in cell_mutations.values():
        for mutation in mutations:
            if mutation not in mutation_counts:
                mutation_counts[mutation] = 0
            mutation_counts[mutation] += 1

    # Calculate VAF for each mutation
    vaf = {mutation: count / num_cells for mutation, count in mutation_counts.items()}
    return vaf

def simulate_bulk_genome(vaf, threshold=0.1):
    """Simulate a bulk genome based on VAF and a threshold."""
    bulk_mutations = [mutation for mutation, freq in vaf.items() if freq > threshold]
    return bulk_mutations

# Collect mutations from all cells
cell_mutations = {}
collect_mutations(tree_root, cell_mutations)

# Calculate VAF for each mutation
vaf = calculate_vaf(cell_mutations, num_cells)

# Simulate the bulk genome (considering only mutations with VAF > 0.1)
bulk_genome = simulate_bulk_genome(vaf, threshold=0.1)

def calculate_bulk_vaf(vaf, num_cells):
    """Calculate and display bulk VAFs for each mutation."""
    bulk_vaf = {mutation: freq for mutation, freq in vaf.items() if freq > 0}
    return bulk_vaf

# Calculate bulk VAF for each mutation
bulk_vaf = calculate_bulk_vaf(vaf, num_cells)

# Display bulk VAFs
print("Bulk genome VAFs (positions and VAF):", bulk_vaf)


# Display results
print("Bulk genome mutations (positions):", bulk_genome)
