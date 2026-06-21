#!/usr/bin/env python3

import argparse
import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# Plotting function
# --------------------------------------------------
def plot_runtime_boxplot(time_files, group_labels, output_folder):
    all_data = []

    for pattern, label in zip(time_files, group_labels):
        files = sorted(glob.glob(pattern))
        if not files:
            raise RuntimeError(f"No files matched pattern: {pattern}")

        for log_file in files:
            df = pd.read_csv(
                log_file,
                sep="\t",
                comment="#",
                names=["cycle", "cycle_time", "cumulative_time"]
            )
            total_runtime = df["cumulative_time"].max()
            all_data.append({
                "group": label,
                "runtime": total_runtime
            })

    df_all = pd.DataFrame(all_data)

    # --------------------------------------------------
    # Prepare data per group
    # --------------------------------------------------
    groups = df_all["group"].unique()
    data_per_group = [
        df_all[df_all["group"] == g]["runtime"].values
        for g in groups
    ]

    # --------------------------------------------------
    # Color palette (print-safe, journal style)
    # --------------------------------------------------
    colors = [
        "#4C72B0",  # muted blue
        "#DD8452",  # muted orange
        "#55A868",  # muted green
        "#C44E52",  # muted red
        "#8172B2",  # muted purple
    ]

    # --------------------------------------------------
    # Plot
    # --------------------------------------------------
    # --------------------------------------------------
    # Plot
    # --------------------------------------------------
    plt.figure(figsize=(6.5, 6))
    
    box = plt.boxplot(
        data_per_group,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=2),
        boxprops=dict(edgecolor="black", linewidth=1.5),
        whiskerprops=dict(color="black", linewidth=1.2),
        capprops=dict(color="black", linewidth=1.2)
    )
    
    # Apply colors to boxes
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
    
    # Overlay points with same color as boxes
    for i, (values, color) in enumerate(zip(data_per_group, colors), start=1):
        jitter = np.random.normal(i, 0.04, size=len(values))
        plt.scatter(
            jitter,
            values,
            s=35,
            facecolors=color,
            edgecolors="black",
            linewidths=0.8,
            zorder=3
        )
    # --------------------------------------------------
    # Aesthetics (journal style)
    # --------------------------------------------------
    plt.xticks(range(1, len(groups) + 1), groups, fontsize=11)
    plt.ylabel("Total simulation runtime (seconds)", fontsize=12)
    plt.title("Runtime distribution across simulation repeats", fontsize=13)

    plt.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.6)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    os.makedirs(output_folder, exist_ok=True)

    pdf_out = os.path.join(output_folder, "simulation_runtime_boxplot.pdf")
    png_out = os.path.join(output_folder, "simulation_runtime_boxplot.png")

    plt.tight_layout()
    plt.savefig(pdf_out)
    plt.savefig(png_out, dpi=300)
    plt.close()

    print(f"Saved plots to:\n  {pdf_out}\n  {png_out}")

# --------------------------------------------------
# Main CLI
# --------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Plot cumulative simulation runtimes as grouped boxplots"
    )

    parser.add_argument(
        "--time_files",
        nargs="+",
        required=True,
        help="Glob patterns for cycle_runtime.log files (one per group)"
    )

    parser.add_argument(
        "--group_labels",
        nargs="+",
        required=True,
        help="Labels for each group (same order as --time_files)"
    )

    parser.add_argument(
        "--output_folder",
        required=True,
        help="Directory to save output plots"
    )

    args = parser.parse_args()

    if len(args.time_files) != len(args.group_labels):
        raise ValueError(
            "--time_files and --group_labels must have the same length"
        )

    plot_runtime_boxplot(
        time_files=args.time_files,
        group_labels=args.group_labels,
        output_folder=args.output_folder
    )

# --------------------------------------------------
# Entry point
# --------------------------------------------------
if __name__ == "__main__":
    main()
