import random
import math


class TreeNode:
    def __init__(self, name, time=0, rate_factor=1.0):
        """
        Represents a node in the coalescent tree.
        :param name: Name of the node (e.g., "cell_0", "anc_1").
        :param time: Time of the node (coalescence time).
        :param rate_factor: Reproductive rate factor for selective sweeps.
        """
        self.name = name
        self.time = time
        self.left = None
        self.right = None
        self.mutations = []
        self.rate_factor = rate_factor  # Default is 1.0 (neutral rate).


class CoalescentTree:
    def __init__(self, num_cells, N, genome_length, mutation_rate):
        """
        Initializes the coalescent tree simulation.
        :param num_cells: Number of cells in the population.
        :param N: Effective population size.
        :param genome_length: Length of the genome.
        :param mutation_rate: Mutation rate per base per generation.
        """
        self.num_cells = num_cells
        self.N = N
        self.genome_length = genome_length
        self.mutation_rate = mutation_rate
        self.tree = None

    def make_coalescence_tree(self):
        """
        Simulates the coalescence tree with selective sweeps.
        Faster-dividing lineages are biased by their rate_factor.
        """
        active_nodes = [TreeNode(name=f"cell_{i}") for i in range(self.num_cells)]
        time = 0

        # Assign faster rates to specific lineages (e.g., every 3rd lineage)
        for i, node in enumerate(active_nodes):
            if i % 3 == 0:  # Example: Every 3rd lineage divides 2x faster
                node.rate_factor = 2.0

        # Coalescence process
        while len(active_nodes) > 1:
            # Total rate adjusted by rate factors
            total_rate = sum(node.rate_factor for node in active_nodes) / (2 * self.N)
            time_to_next = -math.log(random.random()) / total_rate
            time += time_to_next

            # Weighted random selection of two nodes for coalescence
            weights = [node.rate_factor for node in active_nodes]
            selected_nodes = random.choices(active_nodes, weights=weights, k=2)
            active_nodes = [node for node in active_nodes if node not in selected_nodes]

            # Create ancestor node
            ancestor = TreeNode(name=f"anc_{len(active_nodes)}", time=time)
            ancestor.left, ancestor.right = selected_nodes
            ancestor.rate_factor = max(selected_nodes[0].rate_factor, selected_nodes[1].rate_factor)
            active_nodes.append(ancestor)

        self.tree = active_nodes[0]

    def _apply_mutations(self, node, branch_length):
        """
        Applies neutral mutations along a branch.
        :param node: The current node.
        :param branch_length: The branch length to the parent node.
        """
        num_mutations = int(branch_length * self.genome_length * self.mutation_rate)
        new_mutations = [random.randint(0, self.genome_length - 1) for _ in range(num_mutations)]
        node.mutations.extend(new_mutations)

    def assign_mutations(self, node=None):
        """
        Recursively assigns mutations to all branches of the tree.
        :param node: The current node being processed.
        """
        if node is None:
            node = self.tree  # Start from the root

        if node.left and node.right:
            left_branch_length = node.time - node.left.time
            right_branch_length = node.time - node.right.time

            # Apply mutations to left and right branches
            self._apply_mutations(node.left, left_branch_length)
            self._apply_mutations(node.right, right_branch_length)

            # Recursively assign mutations to child nodes
            self.assign_mutations(node.left)
            self.assign_mutations(node.right)

    def print_tree(self, node=None, level=0):
        """
        Prints the tree structure for debugging or visualization.
        :param node: The current node being processed.
        :param level: The depth of the node in the tree.
        """
        if node is None:
            node = self.tree

        print(" " * (level * 4) + f"{node.name} (time: {node.time:.2f}, rate: {node.rate_factor})")
        if node.left:
            self.print_tree(node.left, level + 1)
        if node.right:
            self.print_tree(node.right, level + 1)

    def collect_mutations(self, node=None):
        """
        Collects all mutations from the tree into a dictionary.
        :param node: The current node being processed.
        :return: A dictionary of node names and their mutations.
        """
        if node is None:
            node = self.tree

        mutations = {node.name: node.mutations}
        if node.left:
            mutations.update(self.collect_mutations(node.left))
        if node.right:
            mutations.update(self.collect_mutations(node.right))
        return mutations


# Example Usage
if __name__ == "__main__":
    num_cells = 10  # Number of cells
    N = 1000  # Effective population size
    genome_length = 1_000_000  # Length of the genome (1 Mb)
    mutation_rate = 1e-8  # Mutation rate per base per generation

    # Create the coalescent tree
    coalescent_tree = CoalescentTree(num_cells, N, genome_length, mutation_rate)
    coalescent_tree.make_coalescence_tree()

    # Assign mutations
    coalescent_tree.assign_mutations()

    # Print tree structure
    print("Coalescent Tree:")
    coalescent_tree.print_tree()

    # Collect and print mutations
    mutations = coalescent_tree.collect_mutations()
    print("\nMutations:")
    for node, muts in mutations.items():
        print(f"{node}: {len(muts)} mutations")
