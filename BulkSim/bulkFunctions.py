import pysam
import random
import pickle


def read_fasta(fasta_file):
    """Read the reference genome from a FASTA file."""
    with pysam.FastaFile(fasta_file) as fasta:
        sequences = {ref: fasta.fetch(ref) for ref in fasta.references}
    return sequences


def read_vcf(vcf_file):
    """Read SNPs from a VCF file."""
    snps = {}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            pos = record.pos - 1  # Convert to 0-based indexing
            ref = record.ref
            alt = record.alts
            if len(alt) == 1:  # We only handle biallelic SNPs here
                snps[(record.chrom, pos)] = (ref, alt[0])
    return snps


def apply_snps(reference, snps):
    """Create the modified genome by applying SNPs to the reference genome."""
    modified = {}
    for chrom, seq in reference.items():
        seq_list = list(seq)
        for (snp_chrom, pos), (ref, alt) in snps.items():
            if chrom == snp_chrom and seq_list[pos] == ref:
                # Apply the SNP (replace reference with alternate allele)
                seq_list[pos] = alt
        modified[chrom] = ''.join(seq_list)
    return modified


def generate_bulk_genomes(reference, paternal_snps, maternal_snps, num_copies, coalescent_data_file):
    """Generate multiple copies of paternal and maternal genomes based on VAF distribution."""
    # Create modified paternal and maternal genomes
    paternal_genome = apply_snps(reference, paternal_snps)
    maternal_genome = apply_snps(reference, maternal_snps)

    # Load VAF data from coalescent_data.pkl
    with open(coalescent_data_file, 'rb') as f:
        coalescent_data = pickle.load(f)
    vaf_data = coalescent_data['vaf']

    bulk_paternal = []
    bulk_maternal = []

    for i in range(num_copies):
        # Copy original genomes for each copy
        paternal_copy = paternal_genome.copy()
        maternal_copy = maternal_genome.copy()

        # Apply mutations based on VAF
        for mutation, vaf in vaf_data.items():
            chrom, pos = mutation.split(':')
            pos = int(pos)
            if random.random() <= vaf:
                # Apply mutation for the paternal copy
                if chrom in paternal_copy and len(paternal_copy[chrom]) > pos:
                    paternal_copy[chrom] = paternal_copy[chrom][:pos] + 'N' + paternal_copy[chrom][pos + 1:]

                # Apply mutation for the maternal copy
                if chrom in maternal_copy and len(maternal_copy[chrom]) > pos:
                    maternal_copy[chrom] = maternal_copy[chrom][:pos] + 'N' + maternal_copy[chrom][pos + 1:]

        bulk_paternal.append(paternal_copy)
        bulk_maternal.append(maternal_copy)

    return bulk_paternal, bulk_maternal


def save_fasta(sequences, output_file):
    """Save the sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for chrom, seq in sequences.items():
            f.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


# Example usage:
# Assuming reference_genome, pat_snps, mat_snps, and coalescent_data.pkl are available
reference_genome = read_fasta("reference.fasta")
paternal_snps = read_vcf("pat.vcf")
maternal_snps = read_vcf("mat.vcf")
num_copies = 10

bulk_paternal_genomes, bulk_maternal_genomes = generate_bulk_genomes(
    reference_genome, paternal_snps, maternal_snps, num_copies, 'coalescent_data.pkl')

# Save the bulk genomes
for i, genome in enumerate(bulk_paternal_genomes):
    save_fasta(genome, f"paternal_genome_{i+1}.fasta")

for i, genome in enumerate(bulk_maternal_genomes):
    save_fasta(genome, f"maternal_genome_{i+1}.fasta")
