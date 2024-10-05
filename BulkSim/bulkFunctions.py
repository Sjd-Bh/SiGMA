import random
import pickle
import os
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_vcf(vcf_file):
    snp_positions = {}
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, ref, alt = fields[0], int(fields[1]) - 1, fields[3], fields[4]
            snp_positions[(chrom, pos)] = alt
    return snp_positions

def apply_snps_to_reference(reference_file, snp_positions, output_file):
    records = SeqIO.parse(reference_file, "fasta")
    updated_records = []
    
    for record in records:
        sequence = list(record.seq)
        for (chrom, pos), alt in snp_positions.items():
            if chrom == record.id:
                sequence[pos] = alt
        updated_seq = Seq(''.join(sequence))
        updated_record = SeqRecord(updated_seq, id=record.id, description="with SNPs")
        updated_records.append(updated_record)
    
    SeqIO.write(updated_records, output_file, "fasta")

def amplify_genomes(pat_fasta, mat_fasta, output_fasta, output_vcf, mutations, vaf_info, num_copies=10):
    pat_records = list(SeqIO.parse(pat_fasta, "fasta"))
    mat_records = list(SeqIO.parse(mat_fasta, "fasta"))
    
    amplified_records = []
    vcf_lines = []

    # Write VCF header
    vcf_header = """##fileformat=VCFv4.2
##source=AmplifiedGenomeSimulator
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
    vcf_lines.append(vcf_header)

    mutation_counts = defaultdict(int)

    for i in range(num_copies):
        for record in pat_records + mat_records:
            sequence = list(record.seq)
            chrom = record.id
            for mut in mutations:
                if random.random() < vaf_info.get(mut, 0):  # Use VAF from pkl file
                    original_base = sequence[mut]
                    alt_base = random.choice([b for b in 'ATCG' if b != original_base])
                    sequence[mut] = alt_base
                    mutation_counts[mut] += 1

                    # Write mutation to VCF
                    vcf_line = f"{chrom}\t{mut + 1}\t.\t{original_base}\t{alt_base}\t.\tPASS\t.\n"
                    vcf_lines.append(vcf_line)

            updated_seq = Seq(''.join(sequence))
            amplified_record = SeqRecord(updated_seq, id=f"{record.id}_copy_{i}", description="amplified genome")
            amplified_records.append(amplified_record)
    
    # Save amplified genomes to FASTA
    SeqIO.write(amplified_records, output_fasta, "fasta")

    # Save mutations to VCF
    with open(output_vcf, 'w') as vcf_file:
        vcf_file.writelines(vcf_lines)

    # Print or save observed VAF for verification
    print("Observed Mutation Counts:", dict(mutation_counts))
    print("Expected VAFs:", vaf_info)
    for mut, count in mutation_counts.items():
        observed_vaf = count / (num_copies * len(pat_records + mat_records))
        print(f"Mutation {mut}: Observed VAF = {observed_vaf}, Expected VAF = {vaf_info.get(mut, 0)}")

