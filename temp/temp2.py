
import msprime
import numpy as np
import random

# Parameters
num_cells = 100  # Number of cells to simulate
mutation_rate = 1e-8  # Mutation rate per base per generation
sequence_length = 1e4  # Length of the genome

# Step 1: Generate a coalescent tree
# Using msprime to simulate a tree sequence with cells
tree_sequence = msprime.simulate(sample_size=num_cells, length=sequence_length, mutation_rate=mutation_rate, random_seed=42)

# Step 2: Introduce SNVs
# Extract the mutations from the tree sequence
mutations = list(tree_sequence.mutations())

# Dictionary to track SNV occurrences in cells
snv_dict = {}

# Traverse through the mutations and store information
for mutation in mutations:
    site = mutation.site
    node = mutation.node
    snv_dict[site] = tree_sequence.get_num_samples(node)  # Number of cells with this SNV

# Step 3: Generate bulk genome VAF
bulk_vaf = {}

for site, num_cells_with_snv in snv_dict.items():
    vaf = num_cells_with_snv / num_cells
    bulk_vaf[site] = vaf

# Display the bulk genome VAF
print("Bulk Genome VAF:")
for site, vaf in bulk_vaf.items():
    print(f"SNV at site {site}: VAF = {vaf:.2f}")

# Optional: Further processing or visualization of the bulk VAF
