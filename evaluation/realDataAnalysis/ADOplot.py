import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt

def calculate_dropout_percentage(vcf_file):
    # Read VCF file
    vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
    column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    vcf_data.columns = column_names
    
    # Extract VAF values
    sample_info = vcf_data['SAMPLE'].str.split(':', expand=True)
    ad_dp_values = sample_info[1].str.split(',', expand=True)
    ad_values = ad_dp_values[1].astype(int)
    dp_values = sample_info[2].astype(int)
    vaf_values = ad_values / dp_values
    
    # Calculate dropout percentage
    dropout_percentage = ((vaf_values == 0) | (vaf_values == 1)).mean() * 100
    return dropout_percentage

def main(pta_folder, mda_folder, output_folder):
    # Initialize lists to store data
    pta_percentages = []
    mda_percentages = []

    # Process PTA folder
    for file in os.listdir(pta_folder):
        if file.endswith(".vcf"):
            vcf_file = os.path.join(pta_folder, file)
            percentage = calculate_dropout_percentage(vcf_file)
            pta_percentages.append(percentage)

    # Process MDA folder
    for file in os.listdir(mda_folder):
        if file.endswith(".vcf"):
            vcf_file = os.path.join(mda_folder, file)
            percentage = calculate_dropout_percentage(vcf_file)
            mda_percentages.append(percentage)

    # Plot box plots
    plt.figure(figsize=(8, 6))
    plt.boxplot([pta_percentages, mda_percentages], labels=['PTA', 'MDA'])
    plt.title('Percentage of Positions with VAF=0 or VAF=1')
    plt.ylabel('Percentage')
    plt.xlabel('Group')
    
    # Save the plot
    output_path = os.path.join(output_folder, "dropout_percentage_boxplot.png")
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate dropout percentage and plot box plot')
    parser.add_argument('--pta_folder', type=str, help='Path to the PTA folder containing VCF files')
    parser.add_argument('--mda_folder', type=str, help='Path to the MDA folder containing VCF files')
    parser.add_argument('--output_folder', type=str, default='output', help='Output folder path to save the plot')
    args = parser.parse_args()

    main(args.pta_folder, args.mda_folder, args.output_folder)