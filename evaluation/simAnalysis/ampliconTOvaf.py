import pandas as pd
import os
import argparse
from Bio import SeqIO
import re
###################################################################
# def parse_fasta(filename, ref_fasta):
#     ref_length = len(next(SeqIO.parse(ref_fasta, "fasta")).seq)
#     sequences = []
#     with open(filename, 'r') as f:
#         current_sequence = None
#         for line in f:
#             line = line.strip()
#             if line.startswith('>'):
#                 if current_sequence:
#                     sequences.append(current_sequence)
#                 header = line[1:]
#                 if "Start:" in header and "End:" in header:
#                     start_index = int(header.split("Start: ")[1].split(",")[0])
#                     end_index = int(header.split("End: ")[1])
#                     current_sequence = {'header': header, 'sequence': '', 'start_index': start_index, 'end_index': end_index}
#                 else:
#                     start_index = 1
#                     end_index = ref_length  # Assume end index to be length of reference sequence
#                     current_sequence = {'header': header, 'sequence': '', 'start_index': start_index, 'end_index': end_index}
#             else:
#                 current_sequence['sequence'] += line
#         if current_sequence:
#             sequences.append(current_sequence)
#     return sequences

###################################################################
def parse_fasta(filename):
    sequences_info = []
    for record in SeqIO.parse(filename, "fasta"):
        header = record.description
        start = int(header.split('Start: ')[1].split(',')[0])
        end = int(header.split('End: ')[1])
        sequence = str(record.seq)
        sequences_info.append((header, sequence, start, end))
    return sequences_info

###################################################################
def count_non_reference_snps_df(sequences, snp_positions, ref_nucleotides):
    positions = []
    non_ref_counts = []
    total_counts = []
    vafs = []
    nuc = []
    for position, ref_nucleotide in zip(snp_positions, ref_nucleotides):
        snp_counts = 0
        total_count = 0
        for header, sequence, start, end in sequences:
            if start <= position <= end:
                nucleotide = sequence[position - start]
                print(nucleotide)
                print(header)
                # print(nucleotide)
                nuc.append(nucleotide)
                total_count += 1
                if nucleotide != ref_nucleotide:
                    snp_counts += 1
                
        vaf = snp_counts / total_count if total_count > 4 else 0   
        if vaf < 1:
            positions.append(position)
            non_ref_counts.append(snp_counts)
            total_counts.append(total_count)
            vafs.append(vaf)
            
    snp_counts_df = pd.DataFrame({
        'POS': positions,
        'DP': total_counts,
        'VAF': vafs
    })
    # snp_counts_df = snp_counts_df[(snp_counts_df!=0).all(axis=1)]
    return snp_counts_df

###################################################################
def extract_reference_nucleotides(reference_fasta_file, snp_positions):
    ref_nucleotides = []
    with open(reference_fasta_file, 'r') as f:
        sequence = ''
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line
    for position in snp_positions:
        reference_nucleotide = sequence[position]
        ref_nucleotides.append(reference_nucleotide)
    return ref_nucleotides

###################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process bulk and single-cell data.')
    parser.add_argument('--pat_file', type=str, help='')
    parser.add_argument('--mat_file', type=str, help='')
    parser.add_argument('--single_cell_files', nargs='+', type=str, help='List of paths to single-cell data files')
    parser.add_argument('--output_file', type=str, help='output file path and name')
    parser.add_argument('--ref', type=str, help='reference file')
    parser.add_argument('--chr', type=str, help='output file path and name')
    parser.add_argument('--key1', type=str, help='MDA or COCAllow')
    parser.add_argument('--key2', type=str, help='PTA or COCnotAllow')
    args = parser.parse_args()

    key1 = args.key1
    key2 = args.key2
    patSNP = pd.read_csv(args.pat_file, sep='\t', comment='#', header=None)
    matSNP = pd.read_csv(args.mat_file, sep='\t', comment='#', header=None)

    # merge the dataframe of patSNP and matSNP
    bulk_data = pd.concat([patSNP,matSNP])
    
    column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    bulk_data.columns = column_names
    snp_positions = sorted(bulk_data['POS'])
    merged_data = pd.DataFrame(snp_positions, columns=['POS'])
    single_cell_files = args.single_cell_files
    reference_fasta_file = args.ref
    
    for single_cell_file in single_cell_files:
        sequences = parse_fasta(single_cell_file)
        ref_nucleotides = extract_reference_nucleotides(reference_fasta_file, snp_positions)
        single_cell_data= count_non_reference_snps_df(sequences, snp_positions, ref_nucleotides)
        
        # get the PTA or MDA suffix
        if key1 in single_cell_file:
            prefix = key1
        elif key2 in single_cell_file:
            prefix = key2
        else:
            print(f"prefix {key1} or {key2} is not found in",single_cell_file)
            break
        
      
        # # get the sim suffix
        match_sim = re.search(r'(sim\d+)', single_cell_file)
        if match_sim:
            sim_part = match_sim.group(1)
            
        if 'MDA' in single_cell_file:
            amp = 'MDA'
        elif 'PTA' in single_cell_file:
            amp = 'PTA'

        # # get the colnames
        single_cell_data.columns = [f"{col}_{amp}_{prefix}_{sim_part}" if col != 'POS' else col for col in single_cell_data.columns]
        # single_cell_data.columns = [f"{col}_{prefix}" if col != 'POS' else col for col in single_cell_data.columns]
        
        # Merge single cell data with existing merged data
        merged_data = pd.merge(merged_data, single_cell_data, on='POS', how='left')

        # Replace NaN values in the newly added columns with zeros
        merged_data.dropna(inplace=True)

    # Save the final merged data
    output_file = args.output_file
    merged_data.to_csv(output_file, index=False, sep='\t')
    print(f"Merged data saved to {output_file}")    