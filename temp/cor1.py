import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

# Assuming you have VAF data in a pandas DataFrame where:
# - 'position' column represents SNP positions
# - 'VAF_<sample>' columns represent VAFs for different samples or cells

def calculate_vaf_correlation_by_distance(df, distance_ranges):
    """
    Calculate correlations between SNP VAFs for specified distance ranges.
    
    Parameters:
    df: pandas DataFrame containing SNP positions and VAF values
    distance_ranges: List of distance thresholds to categorize SNP pairs
    
    Returns:
    correlations_by_range: A dictionary where each key is a distance range and
                           each value is a list of correlations for that range
    """
    positions = df['position']
    vaf_cols = df.drop('position', axis=1).columns
    correlations_by_range = {dist: [] for dist in distance_ranges}
    
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            # Calculate genomic distance
            distance = abs(positions[i] - positions[j])
            
            # Calculate correlation between VAFs of the two SNPs
            snp1_vafs = df.iloc[i][vaf_cols]
            snp2_vafs = df.iloc[j][vaf_cols]
            corr, _ = stats.pearsonr(snp1_vafs, snp2_vafs)
            
            # Assign correlation to the appropriate distance range
            for dist in distance_ranges:
                if distance <= dist:
                    correlations_by_range[dist].append(corr)
                    break  # Move to the next SNP pair once the range is matched
    
    return correlations_by_range

def plot_correlations_by_distance(correlations_by_range):
    """
    Plot the distribution of correlations for each distance range.
    
    Parameters:
    correlations_by_range: Dictionary with distance ranges and their correlations
    """
    data = [correlations_by_range[dist] for dist in correlations_by_range]
    labels = [f'<= {dist} bp' for dist in correlations_by_range]
    plt.boxplot(data, labels=labels)
    plt.ylabel('VAF Correlation')
    plt.title('VAF Correlations by Genomic Distance')
    plt.show()

# Load your data (replace with your actual data)
df = pd.DataFrame({
    'position': [100, 200, 350, 800, 1200, 1500, 1800, 2200],
    'VAF_cell_1': [0.45, 0.50, 0.48, 0.15, 0.10, 0.12, 0.20, 0.25],
    'VAF_cell_2': [0.46, 0.52, 0.49, 0.18, 0.12, 0.14, 0.22, 0.27],
    'VAF_cell_3': [0.44, 0.51, 0.47, 0.17, 0.11, 0.13, 0.21, 0.26]
})

# Specify the distance ranges you want to test (e.g., 400, 1000, 1500, 5000 bp)
distance_ranges = [400, 1000, 1500, 5000]

# Calculate correlations
correlations_by_range = calculate_vaf_correlation_by_distance(df, distance_ranges)

# Plot the results
plot_correlations_by_distance(correlations_by_range)

# Optional: Statistical tests to compare correlations between ranges
for dist1, dist2 in zip(distance_ranges[:-1], distance_ranges[1:]):
    corr1 = correlations_by_range[dist1]
    corr2 = correlations_by_range[dist2]
    t_stat, p_value = stats.ttest_ind(corr1, corr2)
    print(f"T-test between {dist1} bp and {dist2} bp: t-stat = {t_stat}, p-value = {p_value}")
