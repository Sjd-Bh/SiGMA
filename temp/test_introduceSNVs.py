from Bio import SeqIO

def get_nucleotide(fasta_file, chrom, position):
    """Retrieve the nucleotide at a given position (0-based index) in a FASTA file."""
    fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    position = position -1 
    if chrom in fasta_dict:
        seq = fasta_dict[chrom]
        if 0 <= position < len(seq):
            return seq[position]
        else:
            return f"Position {position} is out of range for {chrom}."
    else:
        return f"Chromosome {chrom} not found in FASTA file."


def compare_genomes(reference_fasta, mutated_fasta):
    """Compare reference and mutated genome to find actual mutated positions."""
    reference_genome = load_reference_genome(reference_fasta)
    mutated_genome = load_reference_genome(mutated_fasta)

    print("\n### Mutated Positions Report ###")
    for chrom in reference_genome:
        ref_seq = reference_genome[chrom]
        mut_seq = mutated_genome[chrom]
        
        print(f"Checking chromosome: {chrom}")
        for pos in range(len(ref_seq)):
            if ref_seq[pos] != mut_seq[pos]:  # Find positions where bases differ
                print(f"Mutation at position {pos+1}: {ref_seq[pos]} → {mut_seq[pos]}")
    print("\nMutation verification complete.")
