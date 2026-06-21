import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from pathlib import Path
import sys

def setup_cell_press_style():
    """Configures matplotlib for Cell Press publication standards."""
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial']
    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['legend.fontsize'] = 10
    mpl.rcParams['legend.title_fontsize'] = 12
    mpl.rcParams['mathtext.default'] = 'regular' 

def parse_metadata(filename, sample_name):
    # Convert to string to avoid AttributeError on ints/floats
    filename_lower = str(filename).lower()
    sample_lower = str(sample_name).lower()
    
    # Detect Amplification
    if 'pta' in filename_lower or 'pta' in sample_lower:
        amp = 'PTA'
    elif 'mda' in filename_lower or 'mda' in sample_lower:
        amp = 'MDA'
    else:
        amp = 'Unknown'

    # Detect Method
    if 'sccaller' in sample_lower or 'sccaller' in filename_lower:
        method = 'scCaller'
    # Added 'bcftools' to the check below so it gets labeled as BVC
    elif 'bvc' in sample_lower or 'bvc' in filename_lower or 'bcftools' in sample_lower or 'bcftools' in filename_lower:
        method = 'BVC'
    elif 'prosolo' in sample_lower or 'prosolo' in filename_lower:
        method = 'ProSolo'
    else:
        # Fallback to splitting sample name if it's an unrecognized tool
        parts = str(sample_name).split('_')
        method = parts[1] if len(parts) > 1 else 'UnknownMethod'
        
    return method, amp

def main():
    parser = argparse.ArgumentParser(description="Generate publication-quality variant caller plots separated by Amplification.")
    parser.add_argument('--input-files', nargs='+', required=True, help="List of input TSV files")
    parser.add_argument('--output-folder', required=True, help="Directory to save the output PDF plots")
    
    args = parser.parse_args()
    
    out_dir = Path(args.output_folder)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    df_list = []
    for filepath in args.input_files:
        path = Path(filepath)
        if not path.exists():
            continue
            
        temp_df = pd.read_csv(path, sep='\t')
        
        methods, amps = [], []
        for sample in temp_df['Sample']:
            m, a = parse_metadata(path.name, sample)
            methods.append(m)
            amps.append(a)
            
        temp_df['Method'] = methods
        temp_df['Amplification'] = amps
        df_list.append(temp_df)
        
    if not df_list:
        print("No valid data loaded. Exiting.")
        sys.exit(1)
        
    df = pd.concat(df_list, ignore_index=True)
    
    setup_cell_press_style()
    
    # --- Define specific shiny colors for tools ---
    tool_colors = {
        'BVC': '#e41a1c',       # Light Pastel Green
        'ProSolo': '#4daf4a',   # Light Pastel Pink
        'scCaller': '#377eb8',  # Light Pastel Orange
        'combined': '#957DAD',  # ADDED: Color for 'combined' method
        'UnknownMethod': '#808080' # Gray fallback
    }    
    # Generate separate plots for PTA and MDA
    for amp_type in df['Amplification'].unique():
        df_amp = df[df['Amplification'] == amp_type]
        
        if df_amp.empty:
            continue
            
        print(f"Generating plots for {amp_type}...")
        
        # --- Plot 1: Precision vs Recall Scatter ---
        fig, ax1 = plt.subplots(figsize=(5, 5))
        sns.set_style("ticks")
        
        # --- ADDED GRID LINES ---
        ax1.grid(axis='y', linestyle=':', color='lightgray', alpha=0.7, zorder=0)
        
        sns.scatterplot(
            data=df_amp, x='recall', y='precision', 
            hue='Method', palette=tool_colors, 
            s=80, alpha=0.8, zorder=2, ax=ax1  # Added zorder=2 so points sit on top of grid
        )
        
        ax1.set_xlabel('Recall (Sensitivity)')
        ax1.set_ylabel('Precision (PPV)')
        ax1.set_xlim(0, 1.05)
        ax1.set_ylim(0, 1.05)
        ax1.set_title(f'Precision vs Recall ({amp_type})', pad=15, fontweight='bold')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, title='Method')
        sns.despine(ax=ax1)
        
        plt.tight_layout()
        plt.savefig(out_dir / f'Precision_Recall_Scatter_{amp_type}.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        
        # --- Plot 2: F1-Score Boxplot + Stripplot ---
        fig, ax2 = plt.subplots(figsize=(5, 5))
        
        # --- ADDED GRID LINES ---
        ax2.grid(axis='y', linestyle=':', color='lightgray', alpha=0.7, zorder=0)
        
        # Added fill=False so the boxes are transparent but lines match the palette color
        sns.boxplot(
            data=df_amp, x='Method', y='F1score', hue='Method',
            palette=tool_colors, showfliers=False, width=0.5,
            fill=False, legend=False, zorder=1, ax=ax2 
        )
        
        sns.stripplot(
            data=df_amp, x='Method', y='F1score', hue='Method',
            palette=tool_colors, alpha=0.8, linewidth=0, 
            legend=False, zorder=2, ax=ax2
        )
        
        ax2.set_xlabel('Variant Calling Method')
        ax2.set_ylabel('$F_1$ Score') 
        ax2.set_ylim(0, 1.05)
        ax2.set_title(f'$F_1$ Score Comparison ({amp_type})', pad=15, fontweight='bold')
        
        sns.despine(ax=ax2)
        
        plt.tight_layout()
        plt.savefig(out_dir / f'F1_Score_Comparison_{amp_type}.pdf', dpi=300, bbox_inches='tight')
        plt.close()

    print(f"Success! Separate PTA and MDA plots saved to {out_dir.absolute()}")

if __name__ == "__main__":
    main()
