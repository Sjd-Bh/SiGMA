import argparse
from Bio import SeqIO


# Function to mutate the genome sequence at a given position with the alternate base
def mutate_sequence(sequence, position, new_base):
    """
    Replace the base at the specified 1-based position in the sequence with a new base.
    
    Parameters:
    - sequence (str): Original DNA sequence.
    - position (int): 1-based position to mutate.
    - new_base (str): Alternate base to insert.
    
    Returns:
    - str: Mutated DNA sequence.
    """
    pos_index = position - 1  # convert to 0-based index
    sequence_list = list(sequence)
    original_base = sequence_list[pos_index]
    sequence_list[pos_index] = new_base.upper()
    print(f"Mutated position {position}: {original_base} → {new_base}")
    return ''.join(sequence_list)


# Function to read a FASTA file and return the sequence as a string
def read_fasta_sequence(fasta_path):
    """
    Read the first sequence from a FASTA file.
    
    Parameters:
    - fasta_path (str): Path to the reference genome FASTA file.
    
    Returns:
    - str: DNA sequence as a string.
    """
    record = SeqIO.read(fasta_path, "fasta")
    return str(record.seq)


# Function to read the VCF file and return mutations as a dictionary
def read_vcf_file(vcf_path):
    """
    Read a VCF file and return mutations as a dictionary.
    
    Parameters:
    - vcf_path (str): Path to the VCF file.
    
    Returns:
    - dict: {position (int): alt_base (str)}, using 1-based indexing.
    """
    mutations = {}
    with open(vcf_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            position = int(cols[1])  # VCF is 1-based
            alt = cols[4]
            mutations[position] = alt
    return mutations


# Function to apply mutations from the VCF file to the reference genome sequence
def apply_snp_mutations(reference_sequence, vcf_mutations):
    """
    Apply SNP mutations from the VCF file to the reference genome sequence.
    
    Parameters:
    - reference_sequence (str): Original DNA sequence.
    - vcf_mutations (dict): Dictionary of mutations from the VCF file {position: alt_base}.
    
    Returns:
    - str: Mutated genome sequence.
    """
    mutated_sequence = reference_sequence
    for position, alt_base in vcf_mutations.items():
        mutated_sequence = mutate_sequence(mutated_sequence, position, alt_base)
    
    return mutated_sequence


# Function to save the mutated genome to a new FASTA file
def save_mutated_genome(mutated_sequence, output_path):
    """
    Save the mutated genome sequence to a FASTA file.
    
    Parameters:
    - mutated_sequence (str): Mutated genome sequence.
    - output_path (str): Path to save the mutated genome in FASTA format.
    
    Returns:
    - None
    """
    with open(output_path, 'w') as output_file:
        output_file.write(f">mutated_genome\n{mutated_sequence}")
    
    print(f"Mutated genome saved to {output_path}")


# Main function to process the reference genome and apply SNP mutations from the VCF file
def process_genome_with_snp_mutations(fasta_path, vcf_file, output_path):
    """
    Process the reference genome and apply SNP mutations from the VCF file.
    
    Parameters:
    - fasta_path (str): Path to the reference genome FASTA file.
    - vcf_file (str): Path to the VCF file containing SNP mutations.
    - output_path (str): Path to save the mutated genome in FASTA format.
    
    Returns:
    - None
    """
    # Step 1: Read the reference genome sequence
    reference_genome = read_fasta_sequence(fasta_path)
    
    # Step 2: Load VCF mutations
    vcf_mutations = read_vcf_file(vcf_file)
    
    # Step 3: Apply SNP mutations from the VCF file to the reference genome
    mutated_genome = apply_snp_mutations(reference_genome, vcf_mutations)
    
    # Step 4: Save the mutated genome to a new FASTA file
    save_mutated_genome(mutated_genome, output_path)


# Command-line interface using argparse
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Introduce SNP mutations from a VCF file to a reference genome.")
    
    # Define command-line arguments
    parser.add_argument('--ref', required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument('--vcf', required=True, help="Path to the SNP VCF file.")
    parser.add_argument('--out', required=True, help="Path to save the mutated genome in FASTA format.")
    
    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Process the genome and apply SNP mutations
    process_genome_with_snp_mutations(args.ref, args.vcf, args.out)


if __name__ == "__main__":
    main()