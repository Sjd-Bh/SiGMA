import argparse
import pickle
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


# Function to load coalescent data from a pickle file
def load_coalescent_data(pkl_file):
    """
    Load coalescent tree data from a pickle file.
    """
    with open(pkl_file, 'rb') as file:
        coalescent_data = pickle.load(file)
    return coalescent_data


# Function to select mutations from a particular cell (based on its mutations in the coalescent tree)
def select_cell_mutations(node_mutations, cell_id):
    """
    Select mutations from a particular cell (based on its mutations in the coalescent tree).
    """
    if cell_id in node_mutations:
        return node_mutations[cell_id]
    else:
        return {}


# Function to apply mutations from a specific cell to the reference genome sequence
def apply_mutations_from_cell(reference_sequence, node_mutations, cell_id, vcf_mutations):
    """
    Apply mutations from a specific cell to the reference genome sequence, considering
    only those mutations that are present in both the cell_mutations and the VCF file.
    
    Parameters:
    - reference_sequence (str): Original DNA sequence.
    - node_mutations (dict): Dictionary containing mutations for each node/cell.
    - cell_id (str): The cell ID whose mutations will be applied to the genome.
    - vcf_mutations (dict): Dictionary of mutations from the VCF file {position: alt_base}.
    
    Returns:
    - str: Mutated genome sequence.
    """
    # Retrieve mutations for the given cell from coalescent data
    cell_mutations = select_cell_mutations(node_mutations, cell_id)
    
    # Apply mutations sequentially, but only if the mutation exists in the VCF file
    mutated_sequence = reference_sequence
    for position in cell_mutations:
        vcf_position = position + 1  # Convert 0-based position to 1-based
        if vcf_position in vcf_mutations:
            alt_base = vcf_mutations[vcf_position]  # Get the alternate base from VCF
            mutated_sequence = mutate_sequence(mutated_sequence, vcf_position, alt_base)
        else:
            print(f"Warning: Position {vcf_position} from cell {cell_id} not found in VCF mutations.")
    
    return mutated_sequence


# Main function to process the reference genome by applying mutations from a specific cell
def process_genome_with_cell_mutations(fasta_path, pkl_file, snv_file, cell_id, output_path):
    """
    Process the reference genome by applying mutations from a specific cell in the coalescent tree
    and mutations from the provided VCF file.
    
    Parameters:
    - fasta_path (str): Path to the reference genome FASTA file.
    - pkl_file (str): Path to the pickle file containing coalescent tree data.
    - snv_file (str): Path to the SNV VCF file.
    - cell_id (str): The cell ID whose mutations will be applied to the genome.
    - output_path (str): Path to save the mutated genome in FASTA format.
    
    Returns:
    - None
    """
    # Step 1: Read the reference genome sequence
    reference_genome = read_fasta_sequence(fasta_path)
    
    # Step 2: Load coalescent tree data and extract node mutations
    coalescent_data = load_coalescent_data(pkl_file)
    node_mutations = coalescent_data['node_mutations']
    
    # Step 3: Load VCF mutations
    vcf_mutations = read_vcf_file(snv_file)
    
    # Step 4: Apply mutations from the selected cell and VCF mutations
    mutated_genome = apply_mutations_from_cell(reference_genome, node_mutations, cell_id, vcf_mutations)
    
    # Step 5: Save the mutated genome to a new FASTA file
    with open(output_path, 'w') as output_file:
        output_file.write(f">mutated_{cell_id}\n{mutated_genome}")
    
    print(f"Mutated genome for {cell_id} saved to {output_path}")


# Command-line interface using argparse
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Introduce SNVs to a reference genome based on coalescent data.")
    
    # Define command-line arguments
    parser.add_argument('--coal', required=True, help="Path to the coalescent tree pickle file.")
    parser.add_argument('--ref', required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument('--snv', required=True, help="Path to the SNV VCF file.")
    parser.add_argument('--cell', required=True, help="Cell ID (e.g., 'cell_1').")
    parser.add_argument('--out', required=True, help="Path to save the mutated genome in FASTA format.")
    
    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Process the genome with mutations from the selected cell
    process_genome_with_cell_mutations(args.ref, args.coal, args.snv, args.cell, args.out)


if __name__ == "__main__":
    main()