import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# import os

# Function to separate VAF and DP data for MDA and PTA
def MDA_PTA_VAF_DP_separation(dp_vaf_data, key1, key2):
    filtered_cols = {
        'VAF_MDA': dp_vaf_data.filter(like=key1).filter(like='VAF'),
        'VAF_PTA': dp_vaf_data.filter(like=key2).filter(like='VAF'),
        'DP_MDA': dp_vaf_data.filter(like=key1).filter(like='DP'),
        'DP_PTA': dp_vaf_data.filter(like=key2).filter(like='DP')
    }

    vaf_mda, vaf_pta, dp_mda, dp_pta = [filtered_cols[key] for key in filtered_cols]
    
    pos_diff_2_rows = dp_vaf_data['POS'].diff(2).fillna(0)  # 2 rows apart
    
    return vaf_mda, vaf_pta, dp_mda, dp_pta, pos_diff_2_rows

# Function to calculate VAF differences for given distance ranges
# Optimized function to calculate VAF differences
def calculate_vaf_differences_optimized(data, vaf_cols, distance_ranges):
    all_differences = []
    
    # Compute differences between adjacent rows
    data['diff_pos'] = data['POS'].diff().abs()

    for sim in vaf_cols:
        print("loop1")
        for min_dist, max_dist in distance_ranges:
            print("loop2")
            # Filter data within the distance range
            relevant_data = data[(data['diff_pos'] >= min_dist) & (data['diff_pos'] <= max_dist)]

            # Calculate VAF differences with the shift for adjacent rows
            vaf_diff = abs(relevant_data[sim] - relevant_data[sim].shift(-1))

            # Collect valid differences
            valid_differences = vaf_diff.dropna()

            # Append the differences with metadata
            all_differences.append(pd.DataFrame({
                'VAF_Difference': valid_differences,
                'Distance_Range': '{}-{}'.format(min_dist, max_dist),
                'Simulation': sim.split('_')[1],
                'Series': sim.split('_')[2]
            }))
    
    # Concatenate all the results into a single DataFrame
    return pd.concat(all_differences, ignore_index=True)



# Function to calculate correlations based on distance range
def calculate_correlations(pos_diff, data, min_val, max_val):
    indexes = np.where((pos_diff > min_val) & (pos_diff < max_val))[0]
    return [np.corrcoef(data[col][indexes], data[col][indexes - 1])[0, 1] for col in data.columns]

# Function to plot VAF differences and correlations
def cor_plot(csv_file, min_values, max_values, key1, key2, output_folder):
    dp_vaf_data = pd.read_csv(csv_file, sep='\t')
    
    vaf_mda, vaf_pta, dp_mda, dp_pta, pos_diff = MDA_PTA_VAF_DP_separation(dp_vaf_data, key1, key2)

    correlations = {'VAF': {'MDA': [], 'PTA': []}}
    for min_val, max_val in zip(min_values, max_values):
        correlations['VAF']['MDA'].extend(
            [{'cor': cor, 'range': f'{min_val}_{max_val}', 'type': 'MDA'} for cor in calculate_correlations(pos_diff, vaf_mda, min_val, max_val)]
        )
        correlations['VAF']['PTA'].extend(
            [{'cor': cor, 'range': f'{min_val}_{max_val}', 'type': 'PTA'} for cor in calculate_correlations(pos_diff, vaf_pta, min_val, max_val)]
        )

    cor_vaf_df = pd.DataFrame(correlations['VAF']['MDA'] + correlations['VAF']['PTA'])
    
    # VAF Differences
    # distance_ranges = [(0, 400), (400, 1000), (1000, 1500), (1500, 5000)]
    # mda_difference_data = calculate_vaf_differences_optimized(dp_vaf_data, vaf_mda.columns, distance_ranges)
    # pta_difference_data = calculate_vaf_differences_optimized(dp_vaf_data, vaf_pta.columns, distance_ranges)
    
    # difference_data = pd.concat([mda_difference_data, pta_difference_data], ignore_index=True)

    # Plot Correlations
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='type', y='cor', hue='range', data=cor_vaf_df, palette="Set1", width=0.8, whis=1.5, showfliers=False)
    plt.xlabel('Correlations across repeat simulations')
    plt.xticks(rotation=45)
    plt.ylabel('Correlations')
    plt.title('Correlations of VAF based on physical distances')
    plt.tight_layout()
    plt.show()

    # Plot VAF Differences
    # plt.figure(figsize=(12, 8))
    # sns.boxplot(x='Distance_Range', y='VAF_Difference', hue='Simulation', data=difference_data)
    # plt.title('VAF Differences by Distance Range for MDA and PTA Simulations')
    # plt.xlabel('Distance Range (bp)')
    # plt.ylabel('VAF Difference')
    # plt.xticks(rotation=45)
    # plt.tight_layout()
    # plt.show()

# Example usage:
cor_plot('../../test\\MDA_PTA_200kb.csv', [0,200,1000,3000], [200,1000,3000,5000], 'MDA', 'PTA', '../../test\\')
