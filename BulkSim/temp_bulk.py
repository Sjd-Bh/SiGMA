import numpy as np
import random

# Simulated reference genome (short example, replace with real sequence)
reference_genome = "ACGT" * 100  # 400 bp genome (example)

# Number of genome copies in bulk
num_copies = 30

# Number of SNPs
num_snps = 60

# Generate random SNP positions (ensure uniqueness)
snp_positions = random.sample(range(len(reference_genome)), num_snps)

# Generate alternate alleles (A, C, G, T different from reference)
def get_alt_allele(ref):
    bases = {'A', 'C', 'G', 'T'}
    bases.remove(ref)  # Remove reference base
    return random.choice(list(bases))  # Pick a random alt base

# Create a dictionary for SNPs: {position: (ref_base, alt_base)}
snp_dict = {pos: (reference_genome[pos], get_alt_allele(reference_genome[pos])) for pos in snp_positions}

# Generate VAF using binomial distribution (n=30, p=0.5)
vaf_counts = np.random.binomial(n=num_copies, p=0.5, size=num_snps)

# Initialize 30 genome copies as lists for mutability
bulk_genomes = [list(reference_genome) for _ in range(num_copies)]

# Introduce mutations based on VAF counts
for i, (pos, (ref, alt)) in enumerate(snp_dict.items()):
    num_alt = vaf_counts[i]  # Number of copies with the alternate allele
    alt_indices = random.sample(range(num_copies), num_alt)  # Pick strands to mutate
    
    for idx in alt_indices:
        bulk_genomes[idx][pos] = alt  # Replace with alt allele

# Convert genome copies back to strings
bulk_genomes = ["".join(seq) for seq in bulk_genomes]

# Save to FASTA file
def save_fasta(output_file, genome_list):
    """Save all genome copies to a FASTA file."""
    with open(output_file, "w") as f:
        for i, genome in enumerate(genome_list):
            f.write(f">Genome_{i+1}\n")
            f.write(genome + "\n")

# Save the bulk genome
save_fasta("bulk_genome.fasta", bulk_genomes)

print("Bulk genome with heterozygous SNPs saved to bulk_genome.fasta")