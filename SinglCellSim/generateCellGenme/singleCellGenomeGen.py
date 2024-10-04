import random

def generate_single_cell_genome_from_fasta(reference_genome, genome_name="Single_Cell_Genome", 
                                           chromosome="1",num_snps=10, num_snvs=5, num_cnvs=3):
    """
    Generate a synthetic single-cell genome with SNPs, SNVs, and CNVs.
    Args:
        fasta_file (str): The reference genome sequence.
        num_snps (int): Number of SNPs to introduce.
        num_snvs (int): Number of SNVs to introduce.
        num_cnvs (int): Number of CNVs to introduce.
    Returns:
        str: Synthetic single-cell genome in FASTA format.
        list: List of dictionaries containing SNP positions, type, and new nucleotides.
        list: List of dictionaries containing SNV positions, type, and new nucleotides.
        list: List of dictionaries containing CNV positions and length.
        """
    
    def introduce_snps(seq, num_snps,chromosome):        
        snps = []
        if all(base in 'ACGT' for base in seq):
            sequence = seq
        else:
            raise ValueError("Invalid DNA sequence provided.")        
        if len(sequence) < num_snps:
            raise ValueError("The reference genome is too short for the requested number of SNPs.")   
        snp_positions = random.sample(range(len(sequence)), num_snps)
        snp_choices = ['A', 'C', 'G', 'T']        
        for pos in snp_positions:
            old_nucleotide = sequence[pos]
            new_nucleotide = random.choice([x for x in snp_choices if x != old_nucleotide])
            snps.append({'chrom': chromosome, 'position': pos, 'type': f'{old_nucleotide}->{new_nucleotide}', 'new_nucleotide': new_nucleotide})            
            # Convert sequence to a list for manipulation
            sequence_list = list(sequence)
            sequence_list[pos] = new_nucleotide
            sequence = ''.join(sequence_list)   
        snps = [{'position': snp['position'], 'type': snp['type'], 'new_nucleotide': snp['new_nucleotide']} for snp in snps]
        return sequence, snps
    
    def introduce_snvs(sequence, num_snvs,chromosome):
        snvs = []
        snv_positions = random.sample(range(len(sequence)), num_snvs)
        snv_choices = ['A', 'C', 'G', 'T']
        for pos in snv_positions:
            old_nucleotide = sequence[pos]
            new_nucleotide = random.choice([x for x in snv_choices if x != old_nucleotide])
            snvs.append({'chrom': chromosome, 'position': pos, 'type': f'{old_nucleotide}->{new_nucleotide}', 'new_nucleotide': new_nucleotide})
            sequence = sequence[:pos] + new_nucleotide + sequence[pos + 1:]
        return sequence, snvs
    
    def introduce_cnvs(sequence, num_cnvs):
        cnvs = []    
        for _ in range(num_cnvs,chromosome):
            if len(sequence) >= 2:
                cnv_length = random.randint(8000, min(10000, len(sequence)))  # Adjusted to ensure valid range
                cnv_start = random.randint(0, len(sequence) - cnv_length)  # Define the start position of CNV
                cnvs.append({'chrom': chromosome, 'start': cnv_start, 'end': cnv_start + cnv_length})
                cnv_sequence = sequence[cnv_start:cnv_start + cnv_length]
                sequence = sequence[:cnv_start] + cnv_sequence + sequence[cnv_start:]       
        return sequence, cnvs
    # Introduce SNPs, SNVs, and CNVs
    modified_genome, snps = introduce_snps(reference_genome, num_snps)
    modified_genome_permanent ,snvs = introduce_snvs(modified_genome, num_snvs)
    modified_genome, cnvs = introduce_cnvs(modified_genome_permanent, num_cnvs)
    # Format the genome in FASTA format
    fasta_string = modified_genome
    return fasta_string, snps, snvs, cnvs


# Save SNPs, SNVs, and CNVs to separate files
def save_to_vcf(changes, filename, reference_genome):
    with open(filename, "w") as file:
        file.write("##fileformat=VCFv4.2\n")
        file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for change in changes:
            file.write(f"CHROM\t{change['position']+1}\t.\t{reference_genome[change['position']]}\t{change['new_nucleotide']}\t.\tPASS\t.\n")

def save_cnvs_to_bed(cnvs, filename, reference_genome):
    with open(filename, "w") as file:
        for cnv in cnvs:
            file.write(f"CHROM\t{cnv['start']}\t{cnv['end']}\tCNV\n")


def save_fasta(sequence, header, filename):
    with open(filename, "w") as file:
        file.write(f">{header}\n")
        file.write(sequence + "\n")


# Assuming you have paternal_genome and maternal_genome strings

# Save paternal genome as a FASTA file
# save_fasta("Paternal_Genome", paternal_genome, os.path.join(output_dir, 'paternal_genome.fasta'))

# # Save maternal genome as a FASTA file
# save_fasta("Maternal_Genome", maternal_genome, os.path.join(output_dir, 'maternal_genome.fasta'))

# # Concatenate the lines to get the sequence
# sequence = ''.join(line.strip() for line in lines[1:])  # Skipping the first line (header)



        
# # Example Usage
# reference_genome = "ATCGATCGATCG"

# synthetic_genome, snps, snvs, cnvs = generate_single_cell_genome(reference_genome, num_snps=2, num_snvs=1, num_cnvs=1)

# # Save SNPs and SNVs to a VCF-like file
# save_to_vcf(snps, "../test/snps.vcf")


# save_to_vcf(snvs, "../test/snvs.vcf")

# # Save CNVs to a BED file
# save_cnvs_to_bed(cnvs, "../test/cnvs.bed")

# # Save synthetic genome to a FASTA file
# save_fasta(synthetic_genome, "../test/synthetic_genome.fasta")
