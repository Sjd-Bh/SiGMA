import pickle

# Load the .pkl file
with open(pkl_file, "rb") as f:
    data = pickle.load(f)

node_mutations = data["node_mutations"]
tree = data["tree"]

# Get leaf nodes (cells)
leaf_nodes = [node for node in tree.nodes if tree.out_degree(node) == 0]
print("Leaf nodes (cell names):", leaf_nodes)

# Choose one (e.g., cell_0)
cell_id = leaf_nodes[0]
cell_mutations = node_mutations[cell_id]

print(f"Mutations (positions) for {cell_id}:", sorted(cell_mutations))
