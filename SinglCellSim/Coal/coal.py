import msprime
import numpy as np
import random
import tskit

# 1. Simulate the coalescent tree using msprime
def simulate_coalescent_tree(sample_size=10, sequence_length=1e6, recombination_rate=1e-8):
    tree_sequence = msprime.sim_ancestry(samples=sample_size, sequence_length=sequence_length, recombination_rate=recombination_rate)
    return tree_sequence

# 2. Add SNVs along the coalescent tree using msprime's mutation model
def simulate_snvs(tree_sequence, mutation_rate=1e-8):
    mutated_tree = msprime.sim_mutations(tree_sequence, rate=mutation_rate)
    return mutated_tree

# 3. Custom CNV model to introduce duplications/deletions in specific regions
def introduce_cnvs(mutated_tree, cnv_rate=1e-6):
    sequence_length = mutated_tree.sequence_length
    num_cells = mutated_tree.num_samples
    cnv_regions = []

    for i in range(num_cells):
        if random.random() < cnv_rate:  # Introduce a CNV with a certain probability
            start = random.randint(0, int(sequence_length * 0.9))
            length = random.randint(1e4, 1e5)  # CNVs of random length (10kb - 100kb)
            cnv_type = random.choice(['duplication', 'deletion'])
            cnv_regions.append({
                'cell': i,
                'start': start,
                'end': min(start + length, sequence_length),
                'type': cnv_type
            })
    return cnv_regions

# 4. Simulate the effect of CNVs on SNVs
def apply_cnvs_to_snvs(mutated_tree, cnv_regions):
    for variant in mutated_tree.variants():
        position = variant.site.position
        for cnv in cnv_regions:
            if cnv['start'] <= position <= cnv['end']:
                if cnv['type'] == 'deletion':
                    print(f"SNV at position {position} in cell {cnv['cell']} is affected by a deletion.")
                    # Set copy number of the SNV to zero (deleted)
                elif cnv['type'] == 'duplication':
                    print(f"SNV at position {position} in cell {cnv['cell']} is affected by a duplication.")
                    # Increase the copy number of the SNV (duplicated)

# 5. Putting it all together
def simulate_single_cell_genome_with_snvs_and_cnvs(sample_size=10, sequence_length=1e6, recombination_rate=1e-8, mutation_rate=1e-8, cnv_rate=1e-6):
    # Simulate coalescent tree
    tree_sequence = simulate_coalescent_tree(sample_size, sequence_length, recombination_rate)
    
    # Simulate SNVs
    mutated_tree = simulate_snvs(tree_sequence, mutation_rate)
    
    # Introduce CNVs
    cnv_regions = introduce_cnvs(mutated_tree, cnv_rate)
    
    # Apply CNVs to SNVs
    apply_cnvs_to_snvs(mutated_tree, cnv_regions)

    return mutated_tree, cnv_regions

# Example run
if __name__ == "__main__":
    sample_size = 10
    sequence_length = 1e6
    recombination_rate = 1e-8
    mutation_rate = 1e-8
    cnv_rate = 1e-6

    mutated_tree, cnv_regions = simulate_single_cell_genome_with_snvs_and_cnvs(sample_size, sequence_length, recombination_rate, mutation_rate, cnv_rate)
    
    # Output CNV regions
    for cnv in cnv_regions:
        print(f"CNV in cell {cnv['cell']} from {cnv['start']} to {cnv['end']}, type: {cnv['type']}")
