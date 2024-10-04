import random
import numpy as np
from collections import defaultdict

# Step 1: Create a coalescent tree
class TreeNode:
    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.parent = None

    def add_child(self, child):
        self.children.append(child)
        child.parent = self

    def is_leaf(self):
        return len(self.children) == 0

# Generate a simple coalescent tree
def generate_tree(num_cells):
    nodes = [TreeNode(f"Cell_{i}") for i in range(num_cells)]
    while len(nodes) > 1:
        node1 = random.choice(nodes)
        nodes.remove(node1)
        node2 = random.choice(nodes)
        nodes.remove(node2)
        
        parent = TreeNode()
        parent.add_child(node1)
        parent.add_child(node2)
        nodes.append(parent)
    
    return nodes[0]  # Root of the tree

# Step 2: Introduce SNVs in the coalescent tree
def introduce_snvs(tree, snv_probability=0.1):
    snv_dict = defaultdict(list)  # SNV id mapped to nodes where SNV is present

    def traverse_tree(node, current_snvs):
        # Introduce new SNVs with given probability
        if random.random() < snv_probability:
            snv_id = f"SNV_{len(snv_dict) + 1}"
            current_snvs.append(snv_id)
            snv_dict[snv_id].append(node)

        # Pass SNVs to children
        for child in node.children:
            traverse_tree(child, current_snvs.copy())

        # If the node is a leaf, record SNVs for that cell
        if node.is_leaf():
            for snv in current_snvs:
                snv_dict[snv].append(node)

    traverse_tree(tree, [])
    return snv_dict

# Step 3: Generate a bulk genome sequence based on SNV frequencies
def generate_bulk_genome(snv_dict, num_cells):
    # Count only leaf nodes for SNV frequency calculation
    snv_frequency = {snv: len([node for node in nodes if node.is_leaf()]) / num_cells for snv, nodes in snv_dict.items()}
    bulk_genome = []

    for snv, freq in snv_frequency.items():
        if freq >= 0.5:  # Example threshold to consider SNV in bulk genome
            bulk_genome.append(snv)

    return bulk_genome

# Example usage
num_cells = 10
coalescent_tree = generate_tree(num_cells)
snv_dict = introduce_snvs(coalescent_tree, snv_probability=0.2)
bulk_genome = generate_bulk_genome(snv_dict, num_cells)

print("SNV Distribution in Cells:")
for snv, nodes in snv_dict.items():
    print(f"{snv}: {[node.name for node in nodes if node.name is not None]}")

print("\nBulk Genome SNVs:", bulk_genome)
