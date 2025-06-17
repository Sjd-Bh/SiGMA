import argparse
from collections import defaultdict

# Load FASTA file into a dictionary
def load_fasta(fasta_file):
    genome = defaultdict(str)
    with open(fasta_file, 'r') as f:
        chrom = None
        for line in f:
            if line.startswith(">"):
                chrom = line.strip().lstrip(">")
            else:
                genome[chrom] += line.strip()
    return genome

# Load SNVs from VCF file into a list of tuples (chrom, pos, ref, alt)
def load_vcf(vcf_file):
    snvs = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom = parts[0]
                pos = int(parts[1]) - 1  # Convert to 0-based position
                ref = parts[3]
                alt = parts[4]
                snvs.append((chrom, pos, ref, alt))
    return snvs

# Replace main nucleotide with alternate allele from SNVs
def replace_with_alt(genome, snvs):
    for chrom, pos, ref, alt in snvs:
        if chrom in genome:
            # Replace reference nucleotide at the specified position with the alternate
            if genome[chrom][pos] == ref:
                genome[chrom] = genome[chrom][:pos] + alt + genome[chrom][pos + 1:]
    return genome

def save_genome(genome, output_dir):
    """Save the genome to files."""
    for chrom, seq in genome.items():
        output_file = f"{output_dir}/{chrom}_modified.fa"
        with open(output_file, "w") as out:
            out.write(f">{chrom}\n{seq}\n")

def main():
    # Set up command line arguments
    parser = argparse.ArgumentParser(description="Replace reference nucleotides with alternate allele (SNVs) in genome")
    parser.add_argument("--ref", required=True, help="Reference genome (FASTA)")
    parser.add_argument("--snv", required=True, help="VCF file with SNVs")
    parser.add_argument("--output", required=True, help="Output folder")
    
    args = parser.parse_args()

    # Load reference genome
    genome = load_fasta(args.ref)

    # Load SNVs (VCF)
    snvs = load_vcf(args.snv)

    # Replace reference nucleotides with alternate allele (alt) from VCF
    genome = replace_with_alt(genome, snvs)

    # Save the modified genome
    save_genome(genome, args.output)

    print(f"Genome with replaced SNVs has been saved to {args.output}")

if __name__ == "__main__":
    main()