#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
import argparse
import os

# ---------------------------
# Parse command-line arguments
# ---------------------------
parser = argparse.ArgumentParser(description="Summarize ASEReadCounter tables: compute VAF, ADO, LDO per cell.")
parser.add_argument("--input-folder", "-i", required=True,
                    help="Folder containing ASEReadCounter .ase.table files")
parser.add_argument("--output-file", "-o", required=True,
                    help="Output TSV file for summary")
parser.add_argument("--min-dp", "-d", type=int, default=6,
                    help="Minimum DP to consider a site observed (default: 6)")
args = parser.parse_args()

input_folder = args.input_folder
output_file = args.output_file
min_DP = args.min_dp

# ---------------------------
# Gather files
# ---------------------------
files = glob.glob(os.path.join(input_folder, "*.ase.table"))
if len(files) == 0:
    raise ValueError(f"No .ase.table files found in {input_folder}")

# ---------------------------
# Process each file
# ---------------------------
summary = []
all_vafs = []

for file in files:
    df = pd.read_csv(file, sep="\t")
    df["DP"] = df["refCount"] + df["altCount"]

    # Compute VAF
    df["VAF"] = np.where(df["DP"] > 0,
                         df["altCount"] / df["DP"],
                         np.nan)

    n_sites = len(df)

    ldo = (df["DP"] == 0).sum()
    low_conf = ((df["DP"] > 0) & (df["DP"] < min_DP)).sum()
    observed = df[df["DP"] >= min_DP]

    ado = ((observed["VAF"] < 0.1) | (observed["VAF"] > 0.9)).sum()
    balanced = ((observed["VAF"] >= 0.1) & (observed["VAF"] <= 0.9)).sum()

    all_vafs.extend(observed["VAF"].dropna())

    summary.append({
        "sample": os.path.basename(file).replace(".ase.table",""),
        "n_bulk_het": n_sites,
        "n_ldo": ldo,
        "n_lowDP": low_conf,
        "n_observed": len(observed),
        "n_ado": ado,
        "n_balanced": balanced,
        "ADO_rate": ado / max(1, len(observed)),
        "LDO_rate": ldo / max(1, n_sites),
        "Median_VAF": observed["VAF"].median()
    })

# ---------------------------
# Save summary
# ---------------------------
summary_df = pd.DataFrame(summary)
summary_df.to_csv(output_file, sep="\t", index=False)

print(f"Processed {len(files)} files. Summary saved to {output_file}")
print(summary_df)
