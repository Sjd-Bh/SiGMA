import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import re

def extract_genome_size_from_filename(filename):
    match = re.search(r'_(\d+\.?\d*)Mb', filename)
    if match:
        size_mb = float(match.group(1))
        return int(size_mb * 1_000_000)
    else:
        raise ValueError(f"Could not extract genome size from filename: {filename}")

def plot_normalized_ado(csv_files, output_path=None):
    all_data = []

    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            filename = os.path.basename(csv_file)
            genome_size = extract_genome_size_from_filename(filename)

            df['ADO_percent'] = df['ADO'] / genome_size * 100

            # Extract labels
            df['GenomeSize'] = f"{genome_size // 1000} kb"
            df['Method'] = filename.split('_')[-1].replace('.csv', '')  # e.g. "MDA"

            all_data.append(df)

        except Exception as e:
            print(f"Error reading {csv_file}: {e}")

    combined_df = pd.concat(all_data, ignore_index=True)

    plt.figure(figsize=(12, 6))
    sns.boxplot(data=combined_df, x='GenomeSize', y='ADO_percent', hue='Method', palette='Set2')

    plt.title('ADO Rate (% of genome) vs Genome Size')
    plt.xlabel('Genome Size')
    plt.ylabel('ADO Rate (%)')
    plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')

    if output_path:
        plt.savefig(output_path, bbox_inches='tight')
        print(f"Saved plot to {output_path}")
    else:
        plt.tight_layout()
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot ADO rate (%) from ADO summary CSVs, normalized by genome length")
    parser.add_argument('--csv-files', nargs='+', required=True, help='List of ADO summary CSV files')
    parser.add_argument('--output', help='Optional path to save the plot')
    args = parser.parse_args()

    plot_normalized_ado(args.csv_files, args.output)

if __name__ == '__main__':
    main()
