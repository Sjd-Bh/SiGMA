import random
import argparse
import os

def generate_random_sequence(length):
    """Generates a random DNA sequence of a given length."""
    return ''.join(random.choices('ATGC', k=length))

def save_to_fasta(sequence, filename, sequence_name):
    """Saves a DNA sequence to a FASTA file."""
    with open(filename, 'w') as fasta_file:
        fasta_file.write(f'>{sequence_name}\n')
        # Write sequence in lines of 70 characters (as per standard FASTA format)
        for i in range(0, len(sequence), 70):
            fasta_file.write(sequence[i:i+70] + '\n')

def main():
    parser = argparse.ArgumentParser(description="Generate random reference genome sequences of specified lengths.")
    parser.add_argument('-length', type=int, nargs='+', required=True, help="Lengths of sequences to generate in kb (e.g., 200 400 600).")
    parser.add_argument('-o', type=str, required=True, help="Output folder to save the FASTA files.")
    args = parser.parse_args()

    # Ensure the output directory exists
    os.makedirs(args.o, exist_ok=True)

    # Generate and save sequences
    for length_kb in args.length:
        length = length_kb * 1000  # Convert kb to base pairs
        sequence = generate_random_sequence(length)
        filename = os.path.join(args.o, f'reference_sequence_{length_kb}kb.fasta')
        sequence_name = f'{length_kb}kb'
        save_to_fasta(sequence, filename, sequence_name)

    print("Random reference genome sequences generated and saved as FASTA files.")

if __name__ == '__main__':
    main()
