import argparse
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from scipy.stats import pearsonr

def load_data(input_path):
    """
    Loads data. If input_path is a directory, it merges all .tsv/.txt files.
    If input_path is a single file, it loads it directly.
    Assumes columns are sample names and rows are genomic windows.
    """
    if os.path.isdir(input_path):
        files = glob.glob(os.path.join(input_path, "*.*sv")) + glob.glob(os.path.join(input_path, "*.txt"))
        df_list = []
        for f in files:
            name = os.path.splitext(os.path.basename(f))[0]
            # Assuming a single column of counts per file if directory is provided
            tmp_df = pd.read_csv(f, sep='\t', header=None)
            df_list.append(pd.Series(tmp_df.iloc[:, -1].values, name=name))
        df = pd.concat(df_list, axis=1)
    else:
        # Load the merged matrix (ignoring chr, start, end columns if present)
        df = pd.read_csv(input_path, sep='\t')
        # Drop non-numeric genomic coordinate columns if they exist
        cols_to_drop = [c for c in ['chr', 'start', 'end', 'chrom'] if c in df.columns]
        df = df.drop(columns=cols_to_drop)
        
    return df.dropna()

def plot_correlation_bubble(df, output_path, alpha=0.05):
    """Calculates correlations and plots the Cell Reports style bubble plot."""
    n = df.shape[1]
    labels = df.columns
    
    corr_mat = np.zeros((n, n))
    pval_mat = np.zeros((n, n))
    
    # Calculate Pearson correlation
    for i in range(n):
        for j in range(n):
            if i == j:
                corr_mat[i, j], pval_mat[i, j] = 1.0, 0.0
            else:
                r, p = pearsonr(df.iloc[:, i], df.iloc[:, j])
                corr_mat[i, j], pval_mat[i, j] = r, p

    # Setup Figure
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    cmap = plt.get_cmap('RdYlBu_r') # Red to Blue colormap
    norm = Normalize(vmin=-1, vmax=1)
    max_size = 500 # Maximum bubble size

    # Plot lower triangle
    for i in range(n):
        for j in range(n):
            if i > j: # Lower triangle
                r = corr_mat[i, j]
                p = pval_mat[i, j]
                
                # Plot only if statistically significant (p < 0.05)
                if p < alpha:
                    size = abs(r) * max_size
                    color = cmap(norm(r))
                    ax.scatter(j, i, s=size, color=color, edgecolors='black', linewidths=0.5)

    # Formatting
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(labels, rotation=90, fontsize=10)
    ax.set_yticklabels(labels, fontsize=10)
    
    ax.set_xlim(-0.5, n - 1.5)
    ax.set_ylim(n - 0.5, 0.5) # Invert Y axis
    
    for spine in ax.spines.values():
        spine.set_visible(False)
        
    ax.grid(False)
    ax.set_aspect('equal')

    # Colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, shrink=0.7)
    cbar.set_label('Pearson Correlation (r)', rotation=270, labelpad=15, fontsize=12)
    cbar.outline.set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    print(f"Plot saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Generate Amplification Recurrence Bubble Plot")
    parser.add_argument("--input-files", required=True, help="Path to merged matrix file OR folder of input files")
    parser.add_argument("--output-path", required=True, help="Output plot path (e.g., recurrence_plot.pdf)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance threshold (default: 0.05)")
    args = parser.parse_args()

    df = load_data(args.input_files)
    plot_correlation_bubble(df, args.output_path, alpha=args.alpha)

if __name__ == "__main__":
    main()
