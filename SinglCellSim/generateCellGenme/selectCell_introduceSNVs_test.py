import argparse
from Bio import SeqIO

def read_fasta_sequence(fasta_path):
    """
    Read a FASTA file and return the sequence as a string.
    
    Parameters:
    - fasta_path (str): Path to the FASTA file.
    
    Returns:
    - str: DNA sequence as a string.
    """
    record = SeqIO.read(fasta_path, "fasta")
    return str(record.seq)

def compare_genomes(reference_genome, mutated_genome):
    """
    Compare two genome sequences (reference vs mutated).
    This function will report the positions where mutations occurred.
    
    Parameters:
    - reference_genome (str): Reference genome sequence.
    - mutated_genome (str): Mutated genome sequence.
    
    Returns:
    - list: List of mutations, each as a tuple (position, ref_base, mutated_base).
    """
    mutations = []
    if len(reference_genome) != len(mutated_genome):
        raise ValueError("The genomes must have the same length.")
    
    for i in range(len(reference_genome)):
        ref_base = reference_genome[i]
        mutated_base = mutated_genome[i]
        
        # If bases are different, it's a mutation
        if ref_base != mutated_base:
            mutations.append((i + 1, ref_base, mutated_base))  # Position is 1-based
    
    return mutations

def report_comparison(reference_genome, mutated_genome, mutations):
    """
    Report the comparison between the reference and mutated genomes.
    
    Parameters:
    - reference_genome (str): Reference genome sequence.
    - mutated_genome (str): Mutated genome sequence.
    - mutations (list): List of mutations (position, ref_base, mutated_base).
    
    Outputs:
    - None: Prints the comparison report.
    """
    print("Genome Comparison Report")
    print(f"Reference Genome Length: {len(reference_genome)}")
    print(f"Mutated Genome Length: {len(mutated_genome)}")
    
    if len(reference_genome) != len(mutated_genome):
        print("Warning: Genomes have different lengths!")
    
    print(f"\nNumber of mutations: {len(mutations)}\n")
    
    for mutation in mutations:
        position, ref_base, mutated_base = mutation
        print(f"Position {position}: {ref_base} → {mutated_base}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Compare a reference genome and a mutated genome.")
    parser.add_argument('--ref', type=str, required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument('--mutated', type=str, required=True, help="Path to the mutated genome FASTA file.")
    
    args = parser.parse_args()
    
    # Step 1: Read the genomes from the FASTA files
    reference_genome = read_fasta_sequence(args.ref)
    mutated_genome = read_fasta_sequence(args.mutated)
    
    # Step 2: Compare the genomes and get mutations
    mutations = compare_genomes(reference_genome, mutated_genome)
    
    # Step 3: Report the results
    report_comparison(reference_genome, mutated_genome, mutations)

if __name__ == "__main__":
    main()