import pandas as pd
import argparse
import os

def process_ase_table(input_file, phase_file, output_file):
    # Read the GATK ASEReadCounter table
    df = pd.read_csv(input_file, sep='\t')
    
    # Read the phase info extracted via bcftools
    # Columns: contig, position, genotype (e.g., 0|1, 1|0)
    phase_df = pd.read_csv(phase_file, sep='\t', header=None, names=['contig', 'position', 'GT'])
    
    # Merge ASE counts with phase information
    df = pd.merge(df, phase_df, on=['contig', 'position'], how='inner')
    
    # 1. Calculate raw alt VAF
    df['raw_alt_vaf'] = df.apply(
        lambda row: row['altCount'] / row['totalCount'] if row['totalCount'] > 0 else 0, 
        axis=1
    )
    
    # 2. Phase the VAF (Calculate VAF for Haplotype 1)
    # If GT is 1|0, Haplotype 1 is Alt, so VAF_hap1 = raw_alt_vaf
    # If GT is 0|1, Haplotype 1 is Ref, so VAF_hap1 = 1 - raw_alt_vaf
    # Note: Unphased (0/1) will be ignored or left as raw alt VAF depending on your preference.
    def calculate_phased_vaf(row):
        if row['GT'] == '1|0':
            return row['raw_alt_vaf']
        elif row['GT'] == '0|1':
            return 1.0 - row['raw_alt_vaf']
        else:
            return row['raw_alt_vaf'] # Fallback for unphased or missing

    df['phased_raw_vaf'] = df.apply(calculate_phased_vaf, axis=1)
    
    # 3. Whole depth
    df['whole_depth'] = df['rawDepth']
    
    # 4. Filtered phased VAF
    df['filter_vafs'] = df.apply(
        lambda row: row['phased_raw_vaf'] if row['whole_depth'] >= 10 else 0, 
        axis=1
    )
    
    # Rename phased_raw_vaf to raw_vafs to maintain compatibility with downstream scripts
    df = df.rename(columns={'phased_raw_vaf': 'raw_vafs'})
    
    # Select columns
    output_df = df[['contig', 'position', 'GT', 'raw_alt_vaf', 'raw_vafs', 'whole_depth', 'filter_vafs']]

    # Save
    output_df.to_csv(output_file, index=False, sep='\t')
    print(f"File saved successfully to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GATK ASE table to extract Phased VAFs.")
    parser.add_argument("--input", required=True, help="Path to the input ase.table file")
    parser.add_argument("--phase", required=True, help="Path to the phased_snps.txt file")
    parser.add_argument("--output", required=True, help="Path to the output folder")
    
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    
    input_basename = os.path.basename(args.input)
    sample_name = input_basename.replace('.ase.table', '').replace('.table', '').replace('.txt', '')
    
    output_filepath = os.path.join(args.output, f"{sample_name}_phased_vafs.tsv")
    process_ase_table(args.input, args.phase, output_filepath)
