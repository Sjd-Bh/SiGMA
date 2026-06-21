import argparse
import os
import glob
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def get_binned_counts(bam_path, bin_size=10000):
    """Extracts read counts from a BAM file into fixed-size genomic bins."""
    print(f"Extracting bins for: {bam_path} (Bin size: {bin_size/1000}kb)...")
    
    # Check if index exists, warn if not
    if not os.path.exists(bam_path + ".bai") and not os.path.exists(bam_path.replace(".bam", ".bai")):
         print(f"  [!] WARNING: No .bai index found for {bam_path}. Pysam requires indexed BAM files. Please run 'samtools index' first.")
    
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    # Take all chromosomes/contigs present in the BAM file
    chroms = list(bam.references)
    
    if not chroms:
        print(f"  [!] ERROR: Could not find any references in {bam_path}. Check your BAM header.")
        bam.close()
        return pd.DataFrame()
        
    data = []
    for chrom in chroms:
        chrom_len = bam.get_reference_length(chrom)
        for start in range(0, chrom_len, bin_size):
            end = min(start + bin_size, chrom_len)
            # Count reads in this bin
            try:
                count = bam.count(chrom, start, end, read_callback="all")
            except ValueError:
                # Catch errors typically caused by missing BAM indexes
                count = 0 
            # We strip 'chr' just for cleaner labels on the plot, if present
            clean_chrom = chrom.replace('chr', '')
            data.append({'chrom': clean_chrom, 'start': start, 'end': end, 'reads': count})
            
    bam.close()
    return pd.DataFrame(data)

def normalize_and_segment(df):
    """Performs basic normalization and simple rolling-median segmentation."""
    df = df[df['reads'] > 0].copy()
    
    median_reads = df['reads'].median()
    if median_reads == 0 or pd.isna(median_reads):
        median_reads = 1 
        
    df['copy_number'] = (df['reads'] / median_reads) * 2.0
    df['segment'] = df['copy_number'].rolling(window=15, center=True, min_periods=1).median()
    return df

def plot_cnv(df, sample_name="Ampli-1", output_file="cnv_plot.pdf"):
    """Generates a publication-quality CNV plot."""
    plt.rcParams.update({'font.size': 12, 'font.family': 'sans-serif', 'pdf.fonttype': 42, 'axes.linewidth': 1})
    fig, ax = plt.subplots(figsize=(15, 4))
    
    df['genomic_pos'] = range(len(df))
    
    chrom_boundaries = df.groupby('chrom')['genomic_pos'].agg(['min', 'max']).reset_index()
    
    # Safe sorting function for mixed standard and non-standard contigs
    def safe_chrom_sort(c):
        if c.isdigit():
            return int(c)
        elif c.upper() == 'X':
            return 9998
        elif c.upper() == 'Y':
            return 9999
        else:
            # Put other non-standard contigs at the end, but sort them alphabetically amongst themselves
            return 10000

    chrom_boundaries['chrom_num'] = chrom_boundaries['chrom'].apply(safe_chrom_sort)
    # Sort primarily by the numeric mapping, then alphabetically for ties (like custom contigs)
    chrom_boundaries = chrom_boundaries.sort_values(['chrom_num', 'chrom'])
    
    xticks = []
    xticklabels = []
    for i, row in chrom_boundaries.iterrows():
        bg_color = '#e0e0e0' if i % 2 == 0 else '#f2f2f2'
        ax.add_patch(patches.Rectangle((row['min'], -1), row['max'] - row['min'], 10, 
                                       facecolor=bg_color, edgecolor='none', zorder=0))
        xticks.append((row['min'] + row['max']) / 2)
        xticklabels.append(row['chrom'])

    colors = np.where(df['copy_number'] > 2.5, '#e41a1c', 
             np.where(df['copy_number'] < 1.5, '#377eb8', '#808080'))

    ax.scatter(df['genomic_pos'], df['copy_number'], c=colors, s=4, alpha=0.8, zorder=2, edgecolors='none')
    ax.plot(df['genomic_pos'], df['segment'], color='black', linewidth=2, zorder=3)

    ax.set_xlim(0, df['genomic_pos'].max())
    ax.set_ylim(-0.5, 8)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=10, rotation=45 if len(xticks)>25 else 0) # Rotate if many custom contigs
    ax.set_xlabel("Chromosome/Contig", fontweight='bold')
    ax.set_ylabel("Copy Number", fontweight='bold')
    
    ax.text(0.02, 0.85, sample_name, transform=ax.transAxes, fontsize=16, fontweight='bold', color='#1b9e77', zorder=4)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_color('#333333')
    ax.spines['left'].set_color('#333333')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved successfully to {output_file}")
    plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate publication-ready CNV plots from BAM files.")
    parser.add_argument("--input", nargs='+', required=True, help="Path to input BAM file(s). Wildcards supported.")
    parser.add_argument("--output-folder", required=True, help="Directory to save the PDF plots.")
    parser.add_argument("--bin-size", type=int, default=10000, help="Bin size in base pairs (default: 10000 = 10kb)")
    
    args = parser.parse_args()
    os.makedirs(args.output_folder, exist_ok=True)
    
    bam_files = []
    for path in args.input:
        bam_files.extend(glob.glob(path))
        
    if not bam_files:
        print("Error: No BAM files found matching the input criteria.")
        exit(1)
        
    print(f"Found {len(bam_files)} BAM file(s) to process.")
    
    for bam_file in bam_files:
        sample_name = os.path.splitext(os.path.basename(bam_file))[0]
        sample_name = sample_name.replace('.recal.sorted', '') 
        output_pdf = os.path.join(args.output_folder, f"{sample_name}_CNV.pdf")
        
        df_counts = get_binned_counts(bam_file, bin_size=args.bin_size)
        if df_counts.empty:
            print(f"Warning: No valid read counts extracted from {bam_file}. Skipping.")
            continue
            
        df_normalized = normalize_and_segment(df_counts)
        if df_normalized.empty:
            print(f"Warning: After filtering empty bins, no data left for {bam_file}. Skipping.")
            continue
            
        plot_cnv(df_normalized, sample_name=sample_name, output_file=output_pdf)
