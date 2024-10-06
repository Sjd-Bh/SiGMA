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
            node.mutations = node.mutations + new_mutations

        # Choose two nodes to coalesce
        chosen_pair = random.sample(active_nodes, 2)

        # Create a new ancestor node
        new_node = TreeNode(f"anc_{next_available_index}")
        new_node.left, new_node.right = chosen_pair
        new_node.time = chosen_pair[0].time + time_to_next_event

        # Update active nodes
        active_nodes.remove(chosen_pair[0])
        active_nodes.remove(chosen_pair[1])
        active_nodes.append(new_node)

        next_available_index += 1

    # Propagate mutations from the root to all descendant leaves
    propagate_mutations(active_nodes[0], [])

    return active_nodes[0]

def collect_mutations(node, mutations):
    """Collect all unique mutations in the tree."""
    if node is None:
        return
    mutations.update(node.mutations)
    collect_mutations(node.left, mutations)
    collect_mutations(node.right, mutations)

def collect_mutation_counts(node, mutation_counts):
    """Count the occurrences of each mutation in the tree."""
    if node is None:
        return
    for mutation in node.mutations:
        mutation_counts[mutation] += 1
    collect_mutation_counts(node.left, mutation_counts)
    collect_mutation_counts(node.right, mutation_counts)

def collect_node_mutations(node, node_mutations):
    """Collect the mutations for each node in the tree."""
    if node is None:
        return
    node_mutations[node.name] = node.mutations
    collect_node_mutations(node.left, node_mutations)
    collect_node_mutations(node.right, node_mutations)

def collect_leaf_nodes(node, leaf_nodes):
    """Collect all leaf nodes in the tree."""
    if node is None:
        return
    if node.left is None and node.right is None:
        leaf_nodes.append(node)
    else:
        collect_leaf_nodes(node.left, leaf_nodes)
        collect_leaf_nodes(node.right, leaf_nodes)

