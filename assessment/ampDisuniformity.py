#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# 1. Setup Command Line Arguments
parser = argparse.ArgumentParser(description="Calculate and plot Amplification Disuniformity (KL-divergence)")
parser.add_argument("-i", "--input", required=True, help="Input read counts file from samtools bedcov")
parser.add_argument("-o", "--output-folder", default=".", help="Output folder path")
parser.add_argument("-p", "--n-pta", type=int, required=True, help="Number of PTA samples")

# Parse arguments; if none provided, it will print help automatically due to required=True
try:
    args = parser.parse_args()
except SystemExit:
    sys.exit(1)

# Ensure output directory exists
os.makedirs(args.output_folder, exist_ok=True)

# 2. Load the data
# samtools bedcov output is tab-separated without headers
counts = pd.read_csv(args.input, sep='\t', header=None)

# Dynamically name columns based on the input argument
n_cols = counts.shape[1]
n_pta = args.n_pta
n_mda = n_cols - 4 - n_pta  # 3 bed cols + 1 bulk + n_pta

if n_mda < 0:
    print(f"Error: Number of PTA samples ({n_pta}) is too large for the {n_cols} columns in the input file.")
    sys.exit(1)

cell_names = ["bulk"] + [f"PTA_{i}" for i in range(1, n_pta + 1)] + [f"MDA_{i}" for i in range(1, n_mda + 1)]
counts.columns = ["chr", "start", "end"] + cell_names

# 3. Convert counts to probability distributions
probs = counts[cell_names].copy()
probs = probs.div(probs.sum(axis=0), axis=1)

# Replace 0 probabilities with a tiny number to avoid log(0) errors
epsilon = 1e-10
probs = probs.replace(0, epsilon)

# 4. Calculate K-L divergence
kl_results_list = []
Q_bulk = probs["bulk"].values

for cell in cell_names[1:]: # Skip "bulk"
    P_cell = probs[cell].values
    
    # KL-Divergence formula
    kl_div = np.sum(P_cell * np.log(P_cell / Q_bulk))
    
    method_name = "PTA" if "PTA" in cell else "MDA"
    kl_results_list.append({"Cell": cell, "Method": method_name, "KL_Divergence": kl_div})

kl_results = pd.DataFrame(kl_results_list)

# 5. Plotting
# Set general style properties to match Gini plot aesthetics
plt.rcParams.update({
    "font.size": 12,
    "axes.linewidth": 1.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

fig, ax = plt.subplots(figsize=(3.5, 5))

labels = ["PTA", "MDA"]
data = [kl_results[kl_results["Method"] == "PTA"]["KL_Divergence"].values,
        kl_results[kl_results["Method"] == "MDA"]["KL_Divergence"].values]

# Gini script colors
color_map = {
    "PTA": "#5e3c99",  
    "MDA": "#e66101"   
}
colors = [color_map[label] for label in labels]

# Boxplot with empty boxes (facecolor='none') and no outliers
bp = ax.boxplot(data, labels=labels, patch_artist=True, widths=0.3, showfliers=False)

# Customize the colors for boxes, whiskers, caps, and medians
for i, color in enumerate(colors):
    # Box borders
    bp['boxes'][i].set_facecolor('none')
    bp['boxes'][i].set_edgecolor(color)
    bp['boxes'][i].set_linewidth(1)
    
    # Whiskers and caps
    bp['whiskers'][i*2].set(color=color, linewidth=1)
    bp['whiskers'][i*2 + 1].set(color=color, linewidth=1)
    bp['caps'][i*2].set(color=color, linewidth=1)
    bp['caps'][i*2 + 1].set(color=color, linewidth=1)
    
    # Median line
    bp['medians'][i].set(color=color, linewidth=1)

# Add the jittered scatter points
for i, (group_data, color) in enumerate(zip(data, colors)):
    x_positions = np.random.normal(i + 1, 0.05, size=len(group_data))
    ax.scatter(x_positions, group_data, color=color, alpha=0.8, s=5, zorder=3)

# Labels and theme modifications
ax.set_ylim(0, 1) # Force Y-axis to 0 to 1
ax.set_ylabel("Amplification disuniformity\n(KL-divergence to bulk)", 
               fontsize=14, labelpad=10)

# Make X-axis text bold
for label in ax.get_xticklabels():
    label.set_color("black")

ax.tick_params(width=1.5, labelsize=12)

# Remove top and right borders
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Statistical comparison (Wilcoxon rank-sum / Mann-Whitney U test)
stat, p_val = mannwhitneyu(data[0], data[1], alternative='two-sided')

if p_val < 2.2e-16:
    p_label = "p < 2.2e-16"
else:
    p_label = f"p = {p_val:.2e}" if p_val < 0.001 else f"p = {p_val:.3f}"

# Put the text dynamically near the top since the axis is bounded 0 to 1
ax.text(1.5, 0.95, p_label, ha='center', va='bottom', fontsize=12)

# Adjust layout to prevent clipping
plt.tight_layout()

# 6. Save the plot
output_file = os.path.join(args.output_folder, "amplification_disuniformity.pdf")
plt.savefig(output_file, dpi=300, bbox_inches="tight")

print(f"Success! Plot saved to: {output_file}")