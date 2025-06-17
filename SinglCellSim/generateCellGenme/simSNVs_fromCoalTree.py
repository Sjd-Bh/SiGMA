import argparse
import pickle
import random
import os
from Bio import SeqIO

def read_pkl(pkl_file):
    """
    Reads a .pkl file and extracts unique mutation positions in 0-based indexing.
    
    Returns:
        mutations_set: Set of unique mutation positions.
    """
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)  # Load the pkl file

    unified_mutations = set()  # To store unique mutation positions

    # Extract mutation positions from 'mutations' key
    for pos_set in data["mutations"].values():
        unified_mutations.update(pos_set)  # Collect unique positions

    return unified_mutations

def read_fasta(fasta_file):
    """Read the reference genome into a dictionary."""
    ref_genome = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        ref_genome[record.id] = str(record.seq)
    return ref_genome

def get_alt_nucleotide(ref_nuc):
    """Choose a random alternate nucleotide different from the reference."""
    bases = {"A", "T", "C", "G"}
    return random.choice(list(bases - {ref_nuc}))

def write_vcf(vcf_file, snvs, ref_genome):
    """Write SNVs to a VCF file (1-based index)."""
    with open(vcf_file, "w") as vcf:
        # VCF Header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Assuming the reference genome has only one chromosome, get the first chromosome ID
        chrom = list(ref_genome.keys())[0]  # Get the first key (chromosome ID)
        
        for pos in sorted(snvs):  # Sort for consistency
            ref_nuc = ref_genome[chrom][pos]  # Access reference sequence for the chromosome
            alt_nuc = get_alt_nucleotide(ref_nuc)  # Generate alternate nucleotide
            vcf.write(f"{chrom}\t{pos + 1}\t.\t{ref_nuc}\t{alt_nuc}\t.\t.\t.\n")  # Convert to 1-based position

def process_snvs(pkl_file, fasta_file, output_folder):
    """Main function to process SNVs and save VCF files."""
    mutations = read_pkl(pkl_file)
    ref_genome = read_fasta(fasta_file)

    paternal_snvs = set()
    maternal_snvs = set()
    merged_snvs = set()

    # Assign SNVs to paternal or maternal genomes randomly
    for pos in sorted(mutations):
        if random.choice([True, False]):  # Randomly assign to paternal or maternal
            paternal_snvs.add(pos)
        else:
            maternal_snvs.add(pos)
        merged_snvs.add(pos)

    # Save to VCF files (with 1-based positions)
    write_vcf(f"{output_folder}/pat_snvs.vcf", paternal_snvs, ref_genome)
    write_vcf(f"{output_folder}/mat_snvs.vcf", maternal_snvs, ref_genome)
    write_vcf(f"{output_folder}/merged_snvs.vcf", merged_snvs, ref_genome)

def main():
    """Parse command-line arguments and run the SNV simulation."""
    parser = argparse.ArgumentParser(description="Simulate SNVs and generate VCF files.")
    
    # Add command-line arguments
    parser.add_argument('--coal', required=True, help="Path to the coalescent .pkl file.")
    parser.add_argument('--ref', required=True, help="Path to the reference FASTA file.")
    parser.add_argument('--out', required=True, help="Output folder to save the VCF files.")
    
    # Parse arguments
    args = parser.parse_args()

    # Create output folder if it doesn't exist
    os.makedirs(args.out, exist_ok=True)
    
    # Process SNVs using the parsed arguments
    process_snvs(args.coal, args.ref, args.out)

if __name__ == "__main__":
    main()