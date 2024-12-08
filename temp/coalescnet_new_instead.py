import random
import math
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import networkx as nx

class TreeNode:
    def __init__(self, name):
        self.left = None
        self.right = None
        self.time = 0.0
        self.name = name
        self.mutations = []  # List of mutation positions
        self.cnvs = []  # List of CNV events (start, end, type)
        self.selected_mutations = []  # List of mutations under selective sweeps

    def propagate_mutations(self, inherited_mutations, inherited_cnvs, inherited_selected_mutations):
        """
        Propagate mutations, CNVs, and selected mutations down to the child nodes.
        """
        self.mutations = inherited_mutations + self.mutations
        self.cnvs = inherited_cnvs + self.cnvs
        self.selected_mutations = inherited_selected_mutations + self.selected_mutations
        
        if self.left:
            self.left.propagate_mutations(self.mutations, self.cnvs, self.selected_mutations)
        if self.right:
            self.right.propagate_mutations(self.mutations, self.cnvs, self.selected_mutations)

class CoalescentTree:
    def __init__(self, num_cells, N, seed, genome_length, mutation_rate, CNV_rate, CNV_size_mean, CNV_size_std, sweep_strength):
        self.num_cells = num_cells
        self.N = N
        self.seed = seed
        self.genome_length = genome_length
        self.mutation_rate = mutation_rate
        self.CNV_rate = CNV_rate
        self.CNV_size_mean = CNV_size_mean
        self.CNV_size_std = CNV_size_std
        self.sweep_strength = sweep_strength
        self.tree = self.make_coalescence_tree()

    def make_coalescence_tree(self):
        """
        Generate a coalescence tree with mutations, CNVs, and selective sweeps.
        """
        random.seed(self.seed)
        np.random.seed(self.seed)

        active_nodes = [TreeNode(f"cell_{i}") for i in range(self.num_cells)]
        for node in active_nodes:
            node.time = 0

        next_available_index = self.num_cells

        while len(active_nodes) > 1:
            total_rate = len(active_nodes) * (len(active_nodes) - 1) / (2 * self.N)
            time_to_next_event = -math.log(random.random()) / total_rate

            # Apply mutations and CNVs
            for node in active_nodes:
                self._apply_mutations(node, time_to_next_event)

            # Select two nodes to coalesce
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

        # Propagate mutations and CNVs from the root to all descendant leaves
        self._propagate_mutations_and_cnvs(active_nodes[0])

        return active_nodes[0]

    def _apply_mutations(self, node, branch_length):
        """
        Apply point mutations and CNVs to a branch of a given length.
        """
        # Apply point mutations (standard mutation process)
        num_mutations = np.random.poisson(branch_length * self.mutation_rate)
        mutations = [random.randint(0, self.genome_length - 1) for _ in range(num_mutations)]
        node.mutations.extend(mutations)

        # Apply CNVs with some probability based on branch length
        num_cnvs = np.random.poisson(branch_length * self.CNV_rate)
        for _ in range(num_cnvs):
            start_pos = random.randint(0, self.genome_length - 1)
            cnv_size = int(np.random.normal(self.CNV_size_mean, self.CNV_size_std))
            cnv_size = max(1, cnv_size)
            end_pos = min(start_pos + cnv_size, self.genome_length)
            cnv_type = random.choice(['duplication', 'deletion'])

            # Record the CNV event in the node
            node.cnvs.append((start_pos, end_pos, cnv_type))

        # Apply selective sweeps (if any) to this branch
        self._apply_selective_sweep(node, branch_length)

    def _apply_selective_sweep(self, node, branch_length):
        """
        Apply selective sweeps to a branch based on a sweep strength.
        """
        if random.random() < self.sweep_strength:
            sweep_position = random.randint(0, self.genome_length - 1)
            node.selected_mutations.append(sweep_position)

    def _propagate_mutations_and_cnvs(self, root_node):
        """
        Propagate mutations and CNVs from the root node down to all descendant leaves.
        """
        self._propagate_mutations(root_node)

    def _propagate_mutations(self, node):
        """
        Propagate mutations from the root down to all descendant leaves.
        """
        if node.left:
            node.left.propagate_mutations(node.mutations, node.cnvs, node.selected_mutations)
        if node.right:
            node.right.propagate_mutations(node.mutations, node.cnvs, node.selected_mutations)

    def collect_mutations(self, node, mutations):
        """Collect all unique mutations in the tree."""
        if node is None:
            return
        mutations.update(node.mutations)
        self.collect_mutations(node.left, mutations)
        self.collect_mutations(node.right, mutations)

    def collect_mutation_counts(self, node, mutation_counts):
        """Count the occurrences of each mutation in the tree."""
        if node is None:
            return
        for mutation in node.mutations:
            mutation_counts[mutation] += 1
        self.collect_mutation_counts(node.left, mutation_counts)
        self.collect_mutation_counts(node.right, mutation_counts)

    def collect_node_mutations(self, node, node_mutations):
        """Collect the mutations for each node in the tree."""
        if node is None:
            return
        node_mutations[node.name] = node.mutations
        self.collect_node_mutations(node.left, node_mutations)
        self.collect_node_mutations(node.right, node_mutations)

    def collect_leaf_nodes(self, node, leaf_nodes):
        """Collect all leaf nodes in the tree."""
        if node is None:
            return
        if node.left is None and node.right is None:
            leaf_nodes.append(node)
        else:
            self.collect_leaf_nodes(node.left, leaf_nodes)
            self.collect_leaf_nodes(node.right, leaf_nodes)



# def plot_coalescent_tree(tree_root, output_path):
#     """
#     Visualizes the coalescent tree using networkx and matplotlib.

#     Args:
#         tree_root: The root node of the tree.
#         output_path: Path to save the tree visualization.
#     """
#     def add_edges(node, graph):
#         if not node.children:
#             return
#         for child in node.children:
#             graph.add_edge(node.name, child.name)
#             add_edges(child, graph)

#     graph = nx.DiGraph()
#     add_edges(tree_root, graph)

#     pos = nx.spring_layout(graph)  # Layout for tree visualization
#     plt.figure(figsize=(12, 8))
#     nx.draw(
#         graph,
#         pos,
#         with_labels=True,
#         node_size=700,
#         node_color="lightblue",
#         font_size=10,
#         font_weight="bold",
#         arrowsize=15,
#     )
#     plt.title("Coalescent Tree", fontsize=16)
#     plt.savefig(output_path)
#     plt.close()


# Parameters for simulation
# num_cells = 10  # Number of cells (leaf nodes)
# N = 1000  # Effective population size
# genome_length = 10000  # Length of the genome
# mutation_rate = 1e-6  # Mutation rate per site per generation
# CNV_rate = 0.01  # Rate of CNV events
# CNV_size_mean = 1000  # Mean size of CNVs
# CNV_size_std = 500  # Standard deviation of CNV sizes
# sweep_strength = 0.05  # Probability of selective sweep occurring

# # Create coalescent tree with CNVs and selective sweeps
# tree_simulator = CoalescentTree(num_cells, N, seed=42, genome_length=genome_length, mutation_rate=mutation_rate,
#                                 CNV_rate=CNV_rate, CNV_size_mean=CNV_size_mean, CNV_size_std=CNV_size_std,
#                                 sweep_strength=sweep_strength)

# # Collect mutations and CNVs
# mutations = set()
# tree_simulator.collect_mutations(tree_simulator.tree, mutations)
# mutation_counts = defaultdict(int)
# tree_simulator.collect_mutation_counts(tree_simulator.tree, mutation_counts)

# # Collect CNVs and selected mutations
# node_mutations = {} 
# tree_simulator.collect_node_mutations(tree_simulator.tree, node_mutations)
# leaf_nodes = []
# tree_simulator.collect_leaf_nodes(tree_simulator.tree, leaf_nodes)

# # Output the collected mutations and CNVs
# print(f"Total mutations: {len(mutations)}")
# print(f"Mutation counts: {dict(mutation_counts)}")
# print(f"CNVs at leaf nodes: {[(node.name, node.cnvs) for node in leaf_nodes]}")
# print(f"Selected mutations at leaf nodes: {[(node.name, node.selected_mutations) for node in leaf_nodes]}")
