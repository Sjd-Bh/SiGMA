#!/usr/bin/env python3

import argparse
import glob
import os
import numpy as np
import pandas as pd
import configparser
from scipy.stats import pearsonr, spearmanr


# --------------------------------------------------
# Correlation computation
# --------------------------------------------------
def compute_correlations(chrom_data, distance_bins, method):
    # Sort distance bins to use them as upper thresholds
    distance_bins = sorted(distance_bins)
    max_dist = max(distance_bins)
    
    # Dictionaries to hold paired data for each bin
    bins_vaf1 = {b: [] for b in distance_bins}
    bins_vaf2 = {b: [] for b in distance_bins}
    bins_dp1 = {b: [] for b in distance_bins}
    bins_dp2 = {b: [] for b in distance_bins}

    for chrom, records in chrom_data.items():
        records.sort(key=lambda x: x[0])
        n = len(records)
        
        for i in range(n):
            pos1, dp1, vaf1 = records[i]
            # Look ahead for neighboring variants
            for j in range(i + 1, n):
                pos2, dp2, vaf2 = records[j]
                dist = pos2 - pos1
                
                # If distance exceeds our maximum bin, stop searching for pos1
                if dist > max_dist:
                    break
                
                # Find the appropriate bin for this distance
                for b in distance_bins:
                    if dist <= b:
                        bins_vaf1[b].append(vaf1)
                        bins_vaf2[b].append(vaf2)
                        bins_dp1[b].append(dp1)
                        bins_dp2[b].append(dp2)
                        break

    results = {}
    for b in distance_bins:
        # Require at least 30 pairs for a stable correlation calculation
        if len(bins_vaf1[b]) < 30:
            results[b] = (np.nan, np.nan)
            continue

        vaf_1, vaf_2 = bins_vaf1[b], bins_vaf2[b]
        dp_1, dp_2 = bins_dp1[b], bins_dp2[b]

        if method.lower() == "pearson":
            # Log2 normalize depth for Pearson to handle WGA extreme outliers
            dp_1 = np.log2(np.array(dp_1) + 1)
            dp_2 = np.log2(np.array(dp_2) + 1)
            cor_vaf = pearsonr(vaf_1, vaf_2)[0]
            cor_dp = pearsonr(dp_1, dp_2)[0]
        else:
            # Spearman is rank-based, no normalization needed
            cor_vaf = spearmanr(vaf_1, vaf_2)[0]
            cor_dp = spearmanr(dp_1, dp_2)[0]

        results[b] = (cor_vaf, cor_dp)

    return results


# --------------------------------------------------
# Load TSV table
# --------------------------------------------------
def load_tsv(path, min_depth):
    df = pd.read_csv(path, sep="\t")
    chrom_data = {}

    for _, row in df.iterrows():
        dp = row["whole_depth"]
        vaf = row["filter_vafs"]

        # 1. Skip if depth is too low
        # 2. Filter for heterozygous variants ONLY (VAF between 0.1 and 0.9)
        if dp < min_depth or not (0.1 <= vaf <= 0.9):
            continue

        pos = int(row["position"])
        chrom = row["contig"]

        chrom_data.setdefault(chrom, []).append((pos, dp, vaf))

    return chrom_data


# --------------------------------------------------
# Read config.ini
# --------------------------------------------------
def load_config(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)

    min_depth = config.getint("filter", "min_depth")

    distances = [
        int(x.strip())
        for x in config.get("correlation", "distances").split(",")
    ]

    method = config.get("correlation", "method")

    return min_depth, distances, method


# --------------------------------------------------
# Argument parsing
# --------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv-folder", required=True, help="Folder containing the TSV files")
    p.add_argument("--config", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    args = parse_args()

    min_depth, distances, method = load_config(args.config)

    print(f"[CONFIG] min_depth = {min_depth}")
    print(f"[CONFIG] distances (bin thresholds) = {distances}")
    print(f"[CONFIG] method = {method}")

    rows = []

    tsv_files = glob.glob(os.path.join(args.tsv_folder, "*.tsv"))

    for tsv_path in tsv_files:
        cell_id = os.path.basename(tsv_path).replace(".tsv", "").replace("_vafs", "")

        chrom_data = load_tsv(tsv_path, min_depth)
        corrs = compute_correlations(chrom_data, distances, method)

        for d, (cor_vaf, cor_dp) in corrs.items():
            rows.append({
                "cell_id": cell_id,
                "distance_bin_upper_bound": d,
                "cor_vaf": cor_vaf,
                "cor_depth": cor_dp
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output, index=False)

    print("[DONE] Correlation CSV written:", args.output)


if __name__ == "__main__":
    main()
