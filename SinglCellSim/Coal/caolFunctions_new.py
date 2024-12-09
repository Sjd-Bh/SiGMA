import random
import math
import numpy as np
from collections import defaultdict

class TreeNode:
    def __init__(self, name):
        self.left = None
        self.right = None
        self.time = 0.0
        self.name = name
        self.mutations = []

def propagate_mutations(node, inherited_mutations):
    if node is None:
        return
    node.mutations = inherited_mutations + node.mutations
    propagate_mutations(node.left, node.mutations)
    propagate_mutations(node.right, node.mutations)

def finite_site_mutation_model(genome_length, mutation_rate, branch_length):
    """Simulate mutations with a finite site assumption."""
    num_mutations = np.random.poisson(branch_length * mutation_rate * genome_length)
    mutations = [random.randint(0, genome_length - 1) for _ in range(num_mutations)]
    return mutations

def simulate_selective_sweep(node, sweep_probability=0.2):
    """Simulate the effect of a selective sweep by reducing diversity in a subtree."""
    if random.random() < sweep_probability:
        # Remove all mutations from this node and its descendants
        node.mutations = []
    if node.left:
        simulate_selective_sweep(node.left, sweep_probability)
    if node.right:
        simulate_selective_sweep(node.right, sweep_probability)

def make_coalescence_tree(num_cells, N, seed, genome_length, mutation_rate):
    random.seed(seed)
    np.random.seed(seed)

    active_nodes = [TreeNode(f"cell_{i}") for i in range(num_cells)]
    for node in active_nodes:
        node.time = 0

    next_available_index = num_cells

    while len(active_nodes) > 1:
        # Total coalescent rate for the current number of active nodes
        total_rate = len(active_nodes) * (len(active_nodes) - 1) / (2 * N)
        time_to_next_event = np.random.exponential(scale=1 / total_rate)

        # Accumulate mutations along branches
        for node in active_nodes:
            branch_length = 4 * N * mutation_rate * time_to_next_event
            node.mutations += finite_site_mutation_model(genome_length, mutation_rate, branch_length)

        # Choose two nodes to coalesce
        chosen_pair = random.sample(active_nodes, 2)

        # Create a new ancestor node
        new_node = TreeNode(f"anc_{next_available_index}")
        new_node.left, new_node.right = chosen_pair
        new_node.time = max(chosen_pair[0].time, chosen_pair[1].time) + time_to_next_event

        # Update active nodes
        active_nodes.remove(chosen_pair[0])
        active_nodes.remove(chosen_pair[1])
        active_nodes.append(new_node)

        next_available_index += 1

    # Propagate mutations from the root to all descendant leaves
    propagate_mutations(active_nodes[0], [])

    # Simulate a selective sweep in the tree
    simulate_selective_sweep(active_nodes[0])

    return active_nodes[0]

def add_cnvs_and_snvs(node, cnv_probability=0.1, cnv_length_mean=100, genome_length=1e6):
    """Add CNVs and SNVs to the tree."""
    if node is None:
        return
    for i, site in enumerate(node.mutations):
        if random.random() < cnv_probability:
            # Add a CNV starting at this site
            cnv_length = int(np.random.exponential(cnv_length_mean))
            cnv_start = site
            cnv_end = min(site + cnv_length, int(genome_length))
            node.mutations.extend(range(cnv_start, cnv_end))

    add_cnvs_and_snvs(node.left, cnv_probability, cnv_length_mean, genome_length)
    add_cnvs_and_snvs(node.right, cnv_probability, cnv_length_mean, genome_length)

# Collect functions remain the same as your original code

def collect_mutations(node, mutations):
    if node is None:
        return
    mutations.update(node.mutations)
    collect_mutations(node.left, mutations)
    collect_mutations(node.right, mutations)

def collect_mutation_counts(node, mutation_counts):
    if node is None:
        return
    for mutation in node.mutations:
        mutation_counts[mutation] += 1
    collect_mutation_counts(node.left, mutation_counts)
    collect_mutation_counts(node.right, mutation_counts)

def collect_node_mutations(node, node_mutations):
    if node is None:
        return
    node_mutations[node.name] = node.mutations
    collect_node_mutations(node.left, node_mutations)
    collect_node_mutations(node.right, node_mutations)

def collect_leaf_nodes(node, leaf_nodes):
    if node is None:
        return
    if node.left is None and node.right is None:
        leaf_nodes.append(node)
    else:
        collect_leaf_nodes(node.left, leaf_nodes)
        collect_leaf_nodes(node.right, leaf_nodes)
