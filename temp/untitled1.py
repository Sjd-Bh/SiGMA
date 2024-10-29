import pandas as pd
import sys
import os
import argparse

def process_bulk_data(bulk_file,min_threshold,max_threshold):
    bulk_data = pd.read_csv(bulk_file, sep='\t', comment='#', header=None)
    column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    bulk_data.columns = column_names
    sample_info = bulk_data['SAMPLE'].str.split(':', expand=True)
    ad_dp_values = sample_info[1].str.split(',', expand=True)
    ad_values = ad_dp_values[1].astype(int)
    dp_values = sample_info[2].astype(int)
    vaf_values = ad_values / dp_values
    # filtered_positions = sample_info[sample_info[0] == '0/1']
    
   # Filter variants based on VAF values and user-specified thresholds
    filtered_indices = bulk_data[(vaf_values >= min_threshold) & (vaf_values <= max_threshold)].index
    filtered_pos = bulk_data.loc[filtered_indices, 'POS']
    filtered_vaf_values = vaf_values[filtered_indices]
    filtered_dp_values = dp_values[filtered_indices]
    return pd.DataFrame({'POS': filtered_pos, 'VAF_bulk': filtered_vaf_values, 'DP_bulk': filtered_dp_values})

def process_single_cell_data(single_cell_file, filtered_pos):
    single_cell_data = pd.read_csv(single_cell_file, sep='\t', comment='#', header=None)
    column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    single_cell_data.columns = column_names
    sample_info_single_cell = single_cell_data['SAMPLE'].str.split(':', expand=True)
    ad_dp_values_single_cell = sample_info_single_cell[1].str.split(',', expand=True)
    ad_values_single_cell = ad_dp_values_single_cell[1].astype(int)
    dp_values_single_cell = sample_info_single_cell[2].astype(int)
    vaf_values_single_cell = ad_values_single_cell / dp_values_single_cell
    filtered_pos_single_cell = single_cell_data.loc[single_cell_data['POS'].isin(filtered_pos), 'POS']
    filtered_vaf_values_single_cell = vaf_values_single_cell[single_cell_data['POS'].isin(filtered_pos)]
    filtered_dp_values_single_cell = dp_values_single_cell[single_cell_data['POS'].isin(filtered_pos)]
    return pd.DataFrame({'POS': filtered_pos_single_cell, 
                         'VAF': filtered_vaf_values_single_cell, 
                         'DP': filtered_dp_values_single_cell}) 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process bulk and single-cell data.')
    parser.add_argument('--bulk_file', type=str, help='Path to the bulk data file')
    parser.add_argument('--min_threshold', type=float, help='Minimum VAF threshold')
    parser.add_argument('--max_threshold', type=float, help='Maximum VAF threshold')
    parser.add_argument('--output_file', type=str, help='output file path and name')
    parser.add_argument('--single_cell_files', nargs='+', type=str, help='List of paths to single-cell data files')
    args = parser.parse_args()

    bulk_data = process_bulk_data(args.bulk_file, args.min_threshold, args.max_threshold)

    single_cell_files = args.single_cell_files
    merged_data = bulk_data.copy()

    for single_cell_file in single_cell_files:
        single_cell_data = process_single_cell_data(single_cell_file, merged_data['POS'])
        
        # Get the filename for suffix
        file_name = os.path.splitext(os.path.basename(single_cell_file))[0]
        folder_prefix = os.path.basename(os.path.dirname(single_cell_file))
        modified_file_name = f"{folder_prefix}_{file_name}"

        # Add suffix to single cell data columns
        single_cell_data.columns = [f"{col}_{modified_file_name}" if col != 'POS' else col for col in single_cell_data.columns]

        # Merge single cell data with existing merged data
        merged_data = pd.merge(merged_data, single_cell_data, on='POS', how='left')

        # Replace NaN values in the newly added columns with zeros
        for col in single_cell_data.columns[1:]:
            if col.startswith(('VAF_', 'DP_')):
                merged_data[col].fillna(0, inplace=True)

    # Save the final merged data
    output_file = args.output_file
    merged_data.to_csv(output_file, index=False, sep='\t')
    print(f"Merged data saved to {output_file}")    


