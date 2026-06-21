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
    mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Liberation Sans', 'sans-serif']
    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['mathtext.default'] = 'regular'

def parse_metadata(filepath_str, sample_name):
    """Extracts Method and Amplification type."""
    # Convert full path and sample to string and lower case
    path_lower = str(filepath_str).lower()
    sample_lower = str(sample_name).lower()
    
    if 'pta' in path_lower or 'pta' in sample_lower:
        amp = 'PTA'
    elif 'mda' in path_lower or 'mda' in sample_lower:
        amp = 'MDA'
    else:
        amp = 'Unknown'

    # Check for the requested variant callers using the full path
    if 'sccaller' in sample_lower or 'sccaller' in path_lower:
        method = 'scCaller'
    elif 'bvc' in sample_lower or 'bvc' in path_lower:
        method = 'BVC'
    elif 'prosolo' in sample_lower or 'prosolo' in path_lower:
        method = 'ProSolo'
    else:
        method = 'Other'
        
    return method, amp

def main():
    parser = argparse.ArgumentParser(description="Generate publication-quality heatmaps comparing variant callers.")
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
            # Pass the full string of the path, not just path.name
            m, a = parse_metadata(str(path), sample)
            methods.append(m)
            amps.append(a)
        temp_df['Method'] = methods
        temp_df['Amplification'] = amps
        df_list.append(temp_df)
        
    if not df_list:
        print("No valid data loaded. Exiting.")
        sys.exit(1)
        
    df = pd.concat(df_list, ignore_index=True)
    
    # Filter to only keep the methods of interest
    target_methods = ['BVC', 'ProSolo', 'scCaller']
    df = df[df['Method'].isin(target_methods)]

    print("Methods found in data:", df['Method'].unique())
    
    setup_cell_press_style()
    
    # Dictionary mapping column names to display names
    metrics = {
        'recall': 'Recall', 
        'precision': 'Precision', 
        'F1score': '$F_1$ Score'
    }
    
    # Create one plot per metric
    for metric_col, metric_display in metrics.items():
        if metric_col not in df.columns:
            continue
            
        # Group by Method and Amplification to calculate the mean
        agg_df = df.groupby(['Method', 'Amplification'])[metric_col].mean().reset_index()
        
        # Pivot so Methods are rows and Amplifications are columns
        pivot_df = agg_df.pivot(index='Method', columns='Amplification', values=metric_col)
        
        # Reindex to ensure consistent ordering if all methods are present
        ordered_methods = [m for m in target_methods if m in pivot_df.index]
        pivot_df = pivot_df.reindex(ordered_methods)
        
        if pivot_df.empty:
            continue

        fig, ax = plt.subplots(figsize=(4, 3))
        
        # Using RdBu_r to match the Red=High, Blue=Low scale
        sns.heatmap(pivot_df, annot=True, cmap="RdBu_r", vmin=0, vmax=1, 
                    fmt=".3f", cbar_kws={'label': metric_display}, ax=ax,
                    linewidths=0.5, linecolor='white')
        
        ax.set_title(f'Mean {metric_display} Comparison', pad=15)
        ax.set_ylabel('') # Remove 'Method' label 
        ax.set_xlabel('') # Remove 'Amplification' label
        
        plt.xticks(rotation=0)
        plt.yticks(rotation=0)
        
        out_path = out_dir / f'Heatmap_Comparison_{metric_col}.pdf'
        plt.savefig(out_path, dpi=300, bbox_inches='tight')
        plt.close()

    print(f"Success! Comparison heatmaps saved to {out_dir.absolute()}")

if __name__ == "__main__":
    main()
