#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# --------------------------------------------------
# Argument parsing
# --------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot boxplots of correlation vs genomic distance"
    )
    parser.add_argument(
        "--csv",
        required=True,
        help="Input CSV file produced by convert_vcf_cor.py"
    )
    parser.add_argument(
        "--output-folder",
        required=True,
        help="Output folder for plots"
    )
    parser.add_argument(
        "--method",
        required=True,
        choices=["pta", "mda"],
        help="Method of data generation to determine plot colors (pta or mda)"
    )
    return parser.parse_args()


# --------------------------------------------------
# Core plotting function
# --------------------------------------------------
def plot_box(
    df,
    y_col,
    ylabel,
    output_path,
    method
):
    # Ensure distance order
    distances = sorted(df["distance_bin_upper_bound"].unique())
    df["distance_bin_upper_bound"] = pd.Categorical(
        df["distance_bin_upper_bound"],
        categories=distances,
        ordered=True
    )

    # Determine base color and create a gradient palette
    # PTA: #DC143C, MDA: #DAA520
    base_color = "#5e3c99" if method == "pta" else "#e66101"
    
    # Create a palette from dark to light, avoiding pure white at the end
    palette = sns.light_palette(base_color, n_colors=len(distances) + 2, reverse=True)[:-2]

    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(6.5, 4))

    # Plot boxes
    sns.boxplot(
        data=df,
        x="distance_bin_upper_bound",
        y=y_col,
        hue="distance_bin_upper_bound",
        palette=palette,
        dodge=False,
        width=0.6,
        showfliers=False,
        boxprops=dict(facecolor="none", linewidth=1.5),
        ax=ax
    )

    # Update box edge colors and line colors (whiskers, caps, medians)
    for i, patch in enumerate(ax.patches):
        color = palette[i % len(palette)]
        patch.set_edgecolor(color)
        patch.set_facecolor('none')

    # Seaborn adds 6 lines per box (whiskers, caps, median)
    if len(ax.lines) == 6 * len(distances):
        for i in range(len(distances)):
            color = palette[i]
            for j in range(6):
                ax.lines[i * 6 + j].set_color(color)
                ax.lines[i * 6 + j].set_linewidth(1.5)

    # Overlay per-cell points
    sns.stripplot(
        data=df,
        x="distance_bin_upper_bound",
        y=y_col,
        hue="distance_bin_upper_bound",
        palette=palette,
        size=4,
        jitter=True,
        alpha=0.7,
        dodge=False,
        ax=ax
    )

    # Y-axis bounds and styling
    ax.set_ylim(-1, 1)
    
    # --- ADDED: Pale dotted horizontal grid lines behind the data ---
    ax.yaxis.grid(True, linestyle=':', color='lightgray')
    ax.set_axisbelow(True) 
    # ----------------------------------------------------------------
    
    ax.axhline(0, linestyle="--", linewidth=1, color="gray")
    ax.set_xlabel("Genomic distance bin upper bound (bp)", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.tick_params(axis="both", labelsize=10)
    ax.tick_params(axis="x", rotation=90)

    # Clean up legend created by hue
    if ax.legend_ is not None:
        ax.legend_.remove()
        
    sns.despine()
    plt.tight_layout()

    plt.savefig(output_path, format="pdf")
    plt.close()


# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    args = parse_args()

    os.makedirs(args.output_folder, exist_ok=True)

    df = pd.read_csv(args.csv)

    required_cols = {"cell_id", "distance_bin_upper_bound", "cor_vaf", "cor_depth"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")

    # ------------------------
    # VAF correlation plot
    # ------------------------
    plot_box(
        df=df,
        y_col="cor_vaf",
        ylabel="Correlation of VAFs",
        output_path=os.path.join(
            args.output_folder,
            "vaf_correlation_vs_distance.pdf"
        ),
        method=args.method
    )

    # ------------------------
    # Depth correlation plot
    # ------------------------
    plot_box(
        df=df,
        y_col="cor_depth",
        ylabel="Correlation of depth",
        output_path=os.path.join(
            args.output_folder,
            "depth_correlation_vs_distance.pdf"
        ),
        method=args.method
    )

    print("[DONE] Box plots saved to:", args.output_folder)


if __name__ == "__main__":
    main()
