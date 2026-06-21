import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_custom_boxplot(ax, data, position, color, width=0.5):
    """Helper function to plot a transparent boxplot with matching point colors"""
    data = data.dropna()
    if len(data) == 0:
        return
        
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
        plt.setp(bp[element], color=color, linewidth=1)
        
    # Make the box face transparent
    for patch in bp['boxes']:
        patch.set_facecolor('none')
        
    # Add jittered points
    x_jitter = np.random.normal(position, 0.05, size=len(data))
    ax.scatter(x_jitter, data, color=color, alpha=0.7, s=5, zorder=3)

def main():
    parser = argparse.ArgumentParser(description="Plot ADO for PTA and MDA (Simulated vs Real)")
    parser.add_argument('--mda-sim', required=True, help="MDA Simulated summary TSV")
    parser.add_argument('--mda-real', required=True, help="MDA Real summary TSV")
    parser.add_argument('--pta-sim', required=True, help="PTA Simulated summary TSV")
    parser.add_argument('--pta-real', required=True, help="PTA Real summary TSV")
    parser.add_argument('--output', default="ADO_boxplot_publication.pdf", help="Output plot file")
    
    args = parser.parse_args()

    # 1. Load data
    try:
        mda_sim = pd.read_csv(args.mda_sim, sep='\t')
        mda_real = pd.read_csv(args.mda_real, sep='\t')
        pta_sim = pd.read_csv(args.pta_sim, sep='\t')
        pta_real = pd.read_csv(args.pta_real, sep='\t')
    except Exception as e:
        print(f"Error reading files: {e}")
        return

    # 2. Filter real MDA and PTA for samples ending with "chr10.sort"
    mda_real = mda_real[mda_real['sample'].str.contains(r'chr10\.sort', na=False)]
    pta_real = pta_real[pta_real['sample'].str.contains(r'chr10\.sort', na=False)]

    # 3. Extract ADO and convert to percentage
    # Handling potential missing files/columns gracefully
    get_ado = lambda df: df['ADO_rate'] * 100 if 'ADO_rate' in df.columns else pd.Series(dtype=float)
    
    ado_mda_sim = get_ado(mda_sim)
    ado_pta_sim = get_ado(pta_sim)
    ado_mda_real = get_ado(mda_real)
    ado_pta_real = get_ado(pta_real)

    # 4. Publication-Ready Plot Aesthetics
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    
    color_pta = "#5e3c99" 
    color_mda = "#e66101" 
    # Create figure with 2 subplots (Simulated and Real)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5), sharey=True)

    # --- Panel 1: Simulated Data ---
    plot_custom_boxplot(ax1, ado_pta_sim, position=1, color=color_pta)
    plot_custom_boxplot(ax1, ado_mda_sim, position=2, color=color_mda)
    
    ax1.set_title("Simulated Data", fontsize=12)
    ax1.set_xticks([1, 2])
    ax1.set_xticklabels(['PTA', 'MDA'], fontweight='bold')
    ax1.set_ylabel('ADO (%)', fontsize=12)
    ax1.set_xlim(0.5, 2.5)

    # --- Panel 2: Real Data ---
    plot_custom_boxplot(ax2, ado_pta_real, position=1, color=color_pta)
    plot_custom_boxplot(ax2, ado_mda_real, position=2, color=color_mda)

    ax2.set_title("Real Data (chr10.sort)", fontsize=12)
    ax2.set_xticks([1, 2])
    ax2.set_xticklabels(['PTA', 'MDA'])
    ax2.set_xlim(0.5, 2.5)

    # Clean up spines, ticks, and add grid/axis limits
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=11)
        
        # --- NEW CHANGES: Set exact 0-100 limits and add background grid lines ---
        ax.set_ylim(0, 100)
        ax.grid(axis='y', linestyle=':', color='lightgray', alpha=0.7, zorder=0)

    plt.suptitle('Allelic Dropout (ADO) Comparison', fontweight='bold', fontsize=14, y=1.02)
    plt.tight_layout()
    
    # Save high-resolution figures
    plt.savefig(args.output, dpi=600, bbox_inches='tight')
    png_out = args.output.replace('.pdf', '.png')
    plt.savefig(png_out, dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {args.output} and {png_out}")

if __name__ == "__main__":
    main()
