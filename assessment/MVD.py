import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def calculate_mvd(filepath):
    """
    Reads a single cell's TSV file and calculates MVD from the filter_vafs column.
    MVD = median(|VAF - 0.5|)
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        valid_vafs = df.loc[df['filter_vafs'] > 0.0, 'filter_vafs'].dropna()
        if len(valid_vafs) == 0:
            return np.nan
        return np.median(np.abs(valid_vafs - 0.5))
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return np.nan

# Changed width from 0.25 to 0.5 to make boxes wider
def plot_custom_boxplot(ax, data, position, color, width=0.5):
    """Helper function to plot a transparent boxplot with matching point colors"""
    # Create the boxplot
    bp = ax.boxplot(
        data, 
        positions=[position], 
        widths=width, 
        patch_artist=True, 
        showfliers=False
    )
    
    # Set colors for all boxplot elements
    for element in ['boxes', 'whiskers', 'caps', 'medians']:
        plt.setp(bp[element], color=color, linewidth=1, zorder=3)
        
    # Make the box face transparent
    for patch in bp['boxes']:
        patch.set_facecolor('none')
        
    # Add jittered points
    x_jitter = np.random.normal(position, 0.04, size=len(data))
    ax.scatter(x_jitter, data, color=color, alpha=0.7, s=5, zorder=4)

def main():
    parser = argparse.ArgumentParser(description="Plot Median VAF Deviance (MVD) and paired differences")
    parser.add_argument('--input-pta', nargs='+', required=True, help="List of PTA TSV files")
    parser.add_argument('--input-mda', nargs='+', required=True, help="List of MDA TSV files")
    parser.add_argument('--output', required=True, help="Output image file (e.g., mvd_plot.pdf)")
    
    args = parser.parse_args()

    # Ensure paired inputs have the same length
    if len(args.input_pta) != len(args.input_mda):
        print("Error: The number of PTA and MDA input files must be exactly the same for paired analysis.")
        return

    # Calculate MVD (assuming pairs are passed in the exact same order)
    pta_mvds = []
    mda_mvds = []
    
    for f_pta, f_mda in zip(args.input_pta, args.input_mda):
        pta_mvd = calculate_mvd(f_pta)
        mda_mvd = calculate_mvd(f_mda)
        if not np.isnan(pta_mvd) and not np.isnan(mda_mvd):
            pta_mvds.append(pta_mvd)
            mda_mvds.append(mda_mvd)
            
    if not pta_mvds:
        print("No valid paired data found to plot.")
        return

    # Calculate paired differences (MDA - PTA)
    mvd_diffs = np.array(mda_mvds) - np.array(pta_mvds)

    # -----------------------------------------
    # Publication-Ready Plot Aesthetics
    # -----------------------------------------
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    
    # Sharp colors: Teal (PTA), Yellow Ocher (MDA), Crimson (Diff)
    color_pta = '#542788'
    color_mda = '#e66101' 
    color_diff = '#b35806'

    # Create figure with 1 plot area
    fig, ax1 = plt.subplots(figsize=(5.5, 4.5))
    
    # Create the secondary Y-axis (Right)
    ax2 = ax1.twinx()

    # Add light background dotted lines on the y-axis
    ax1.grid(axis='y', linestyle=':', color='lightgray', linewidth=1, zorder=0)
    ax1.set_axisbelow(True) # Ensure grid is behind the plot elements

    # --- Plotting MVD Boxplots on Left Axis ---
    plot_custom_boxplot(ax1, pta_mvds, position=1, color=color_pta)
    plot_custom_boxplot(ax1, mda_mvds, position=2, color=color_mda)

    # --- Plotting MVD-diff Boxplot on Right Axis ---
    plot_custom_boxplot(ax2, mvd_diffs, position=3, color=color_diff)
    
    # Add a dashed line at 0 to indicate no difference
    ax2.axhline(0, color='gray', linestyle='--', linewidth=1, zorder=1)

    # Set Axes Labels and Ticks
    ax1.set_xticks([1, 2, 3])
    ax1.set_xticklabels(['PTA', 'MDA', 'MDA - PTA'])
    ax1.set_xlim(0.5, 3.5)
    
    ax1.set_ylabel('Median VAF Deviance (MVD)', fontsize=12)
    ax2.set_ylabel('MVD Difference', fontsize=12)

    # Force identical scales to 0.5 (Starting at -0.1 to leave room for negative differences)
    ax1.set_ylim(-0.1, 0.5)
    ax2.set_ylim(-0.1, 0.5)

    # Clean up top spine for both axes (Right spine must remain visible for ax2)
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax2.tick_params(axis='both', which='major', labelsize=10)

    plt.suptitle('Allelic Distortion & Paired Difference', fontweight='bold', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {args.output}")

if __name__ == "__main__":
    main()
