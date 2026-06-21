import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def plot_publication_vaf(input_files, output_prefix):
    # 1. Setup Publication-Quality Aesthetics
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'axes.linewidth': 1.2,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'legend.title_fontsize': 10
    })

    data_frames = []
    
    # Auto-detect group based on output prefix to assign the correct deep color
    prefix_upper = output_prefix.upper()
    if "PTA" in prefix_upper:
        main_color = "#5e3c99"  # Crimson
        group_label = "PTA"
    else:
        main_color = "#e66101"  # Goldenrod (MDA default)
        group_label = "MDA"

    # 2. Process Input Files
    for file in input_files:
        basename = os.path.basename(file)
        sample_name = basename.replace('_vafs.tsv', '').replace('.tsv', '')
        
        df = pd.read_csv(file, sep='\t')
        
        valid_vafs = df[df['filter_vafs'] > 0].copy()
        valid_vafs['Sample'] = sample_name
        data_frames.append(valid_vafs)
        
    if not data_frames:
        print("Error: No valid data found after filtering.")
        return

    combined_df = pd.concat(data_frames, ignore_index=True)
    
    # 3. Plotting
    # Slightly wider figure to make room for the legend on the right
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.set_style("ticks")
    
    # Generate a rainbow-like palette for the individual samples
    sample_colors = sns.color_palette("hls", len(data_frames))
    
    # Plot individual samples (Pale dotted lines/fills with distinct colors)
    for df, color in zip(data_frames, sample_colors):
        sample_name = df['Sample'].iloc[0]
        sns.kdeplot(
            data=df, 
            x='filter_vafs', 
            fill=True, 
            alpha=0.15,        # Transparency creates the pale shadow effect
            linewidth=1.5, 
            linestyle=':',     # Dotted line
            clip=(0, 1),
            color=color,       # Unique color from the palette
            label=sample_name,
            ax=ax,
            warn_singular=False
        )
    
    # Plot the Main Distribution (Thick Mean Line representing combined data)
    sns.kdeplot(
        data=combined_df, 
        x='filter_vafs', 
        fill=False, 
        linewidth=4,      # Deep, thick line
        clip=(0, 1),
        color=main_color,
        label=f"Main Average ({group_label})",
        ax=ax,
        zorder=10           # Forces the main line to be drawn on top of the fills
    )
    
    # 4. Formatting Axes and Labels
    ax.set_xlabel('Variant Allele Frequency (VAF)')
    ax.set_ylabel('Density')
    
    # Fix X and Y limits
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 5.0)     # Fixed y-axis to 5
    
    ax.set_xticks([0, 0.25, 0.50, 0.75, 1.0])
    
    sns.despine(ax=ax)
    
    # Move the legend outside the plot area
    ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1.05, 0.5), title="Samples")
    
    # 5. Save the Plot
    # bbox_inches='tight' is crucial here so the external legend isn't cut off when saving
    pdf_out = f"{output_prefix}.pdf"
    plt.savefig(pdf_out, format='pdf', dpi=300, bbox_inches='tight')
    
    png_out = f"{output_prefix}.png"
    plt.savefig(png_out, format='png', dpi=600, bbox_inches='tight')
    
    print(f"Publication-ready plots saved as {pdf_out} and {png_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate publication-quality VAF distribution plot.")
    parser.add_argument("--inputs", nargs='+', required=True, help="List of input TSV files (e.g., *_vafs.tsv)")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files (e.g., /path/to/MDA_VAF_Plot)")
    
    args = parser.parse_args()
    plot_publication_vaf(args.inputs, args.output_prefix)
