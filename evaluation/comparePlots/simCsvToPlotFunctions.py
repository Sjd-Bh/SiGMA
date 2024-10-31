###############################################################################
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re

###############################################################################
## separation csv data cols for VAF and depth of PTA and MDA
def MDA_PTA_VAF_DP_separation(dp_vaf_data, key1, key2):
    # Filter columns based on conditions and extract their names
    filtered_cols = {
        'VAF_MDA': dp_vaf_data.filter(like=key1).filter(like='VAF'),
        'VAF_PTA': dp_vaf_data.filter(like=key2).filter(like='VAF'),
        'DP_MDA': dp_vaf_data.filter(like=key1).filter(like='DP'),
        'DP_PTA': dp_vaf_data.filter(like=key2).filter(like='DP')
    }

    # Extract data based on filtered columns
    vaf_mda, vaf_pta, dp_mda, dp_pta = [filtered_cols[key] for key in filtered_cols]
    
    # Calculate POS column differences
    # pos_diff = dp_vaf_data['POS'].diff().fillna(0)
    
    # Calculate POS column differences for 1 row apart and 2 rows apart
    pos_diff_1_row = dp_vaf_data['POS'].diff(1).fillna(0)  # 1 row apart
    pos_diff_2_rows = dp_vaf_data['POS'].diff(2).fillna(0)  # 2 rows apart
    
    # Return the four sets of data and pos_diff as a tuple
    return vaf_mda, vaf_pta, dp_mda, dp_pta, pos_diff_2_rows

###############################################################################
## calculate correlation VAF based on physical distance of SNV sites
def calculate_correlations(pos_diff, data, min_val, max_val):
    indexes = np.where((pos_diff > min_val) & (pos_diff < max_val))[0]
    return [np.corrcoef(data[col][indexes], data[col][indexes - 1])[0, 1] for col in data.columns]

## plot the calculated correlatios (each box is for each min and max distance)
def cor_plot(pat_snp_file,csv_file, min_values, max_values, key1, key2, output_folder):
    # Load data
    dp_vaf_data = pd.read_csv(csv_file, sep='\t')
    filename_prefix = os.path.splitext(os.path.basename(csv_file))[0]

    # Load SNP positions and phase data
    pat_snp = pd.read_csv(pat_snp_file, sep='\t', comment='#', header=None,
                           names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    pat_snp = pat_snp.drop_duplicates(subset=['POS'])
    snp_positions = set(pat_snp['POS'])
    vaf_cols = [col for col in dp_vaf_data.columns if 'VAF' in col]

    # Phase the data
    mask = dp_vaf_data['POS'].isin(snp_positions)
    dp_vaf_data.loc[mask, vaf_cols] = 1 - dp_vaf_data.loc[mask, vaf_cols]

    # Separate data into MDA/PTA for VAF and DP
    vaf_mda, vaf_pta, dp_mda, dp_pta, pos_diff = MDA_PTA_VAF_DP_separation(dp_vaf_data, key1, key2)

    # Calculate correlations
    correlations = {'VAF': {'MDA': [], 'PTA': []}, 'DP': {'MDA': [], 'PTA': []}}
    for min_val, max_val in zip(min_values, max_values):
        correlations['VAF']['MDA'].extend(
            [{'cor': cor, 'range': f'{min_val}_{max_val}', 'type': 'MDA'} for cor in calculate_correlations(pos_diff, vaf_mda, min_val, max_val)]
        )
        correlations['VAF']['PTA'].extend(
            [{'cor': cor, 'range': f'{min_val}_{max_val}', 'type': 'PTA'} for cor in calculate_correlations(pos_diff, vaf_pta, min_val, max_val)]
        )
        correlations['DP']['MDA'].extend(
            [{'cor': cor, 'range': f'{min_val}_{max_val}', 'type': 'MDA'} for cor in calculate_correlations(pos_diff, dp_mda, min_val, max_val)]
        )
        correlations['DP']['PTA'].extend(
            [{'cor': cor, 'range': f'{min_val}_{max_val}', 'type': 'PTA'} for cor in calculate_correlations(pos_diff, dp_pta, min_val, max_val)]
        )

    # Prepare DataFrames for plotting
    cor_vaf_df = pd.DataFrame(correlations['VAF']['MDA'] + correlations['VAF']['PTA'])
    # cor_dp_df = pd.DataFrame(correlations['DP']['MDA'] + correlations['DP']['PTA'])

    # Plot VAF correlations
    output_plot_path = os.path.join(output_folder, f'correlations_DP_{filename_prefix}.png')
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='type', y='cor', hue='range', data=cor_vaf_df, palette="Set1", width=0.8, whis=1.5, showfliers=False)
    plt.xlabel('Correlations across repeat simulations')
    plt.xticks(rotation=45)
    plt.ylabel('Correlations')
    plt.title('Correlations of Depth of coverage based on physical distances')
    plt.savefig(output_plot_path)  # Save the plot
    plt.show()
    print('Boxplot was saved to:', output_plot_path)
###############################################################################
## calculating MVD-Diff (Median VAF Deviance Difference between MDA and PTA) statistic
def calculate_mvd(row, key1, key2):
    vaf_allow = np.abs(row.filter(like=f'VAF_{key1}').values - 0.5)
    vaf_not_allow = np.abs(row.filter(like=f'VAF_{key2}').values - 0.5)
    return pd.Series({'Statistic1': np.median(vaf_allow), 'Statistic2': np.median(vaf_not_allow)})

## Plot MVD_Diff
def mvd_diff_mp(csv_files, key1, key2, output_folder):
    merged_data = []

    for csv_file in csv_files:
        dp_vaf_data = pd.read_csv(csv_file, sep='\t')
        filename_prefix = os.path.splitext(os.path.basename(csv_file))[0]
        
        # Define key pairs to calculate differences
        keys = [(key1, key2)]
        for key_pair in keys:
            # Calculate statistics for the current key pair
            statistics_df = dp_vaf_data.apply(calculate_mvd, args=key_pair, axis=1).dropna()
            statistics_df['Difference'] = statistics_df['Statistic1'] - statistics_df['Statistic2']
            statistics_df['statType'] = f'{key_pair[0]} vs {key_pair[1]}'
            statistics_df['amp'] = filename_prefix

            # Append the data for plotting
            merged_data.append(statistics_df[['Difference', 'statType', 'amp']])

    # Concatenate all data for plotting
    merged_data = pd.concat(merged_data, ignore_index=True)

    # Save the merged data to a CSV file
    output_csv = os.path.join(output_folder, f'wilcoxon_{filename_prefix}.csv')
    merged_data.to_csv(output_csv, index=False, sep='\t')

    # Plot the boxplot of differences for all key pairs
    plt.figure(figsize=(10, 8))
    sns.set_theme(style="ticks", palette="pastel")
    ax = sns.boxplot(x='statType', y='Difference', hue='amp', data=merged_data, palette="Set1", width=0.6, whis=1.5, showfliers=False)
    
    # Calculate the mean values and annotate them above each box
    mean_values = merged_data.groupby(['statType', 'amp'])['Difference'].mean().reset_index()
    for i, row in mean_values.iterrows():
        stat_type = row['statType']
        amp = row['amp']
        mean = row['Difference']
        
        # Find position for the x-tick corresponding to the `statType` and `amp`
        pos = merged_data[(merged_data['statType'] == stat_type) & (merged_data['amp'] == amp)].index[0]
        ax.text(pos, mean + 0.02, f'{mean:.2f}', ha='center', va='bottom', fontsize=9, color='black')

    plt.xlabel('Reference length')
    plt.xticks(rotation=45)
    plt.ylabel('MVD-Diff')
    plt.title('Median VAF Deviance Difference between MDA and PTA')
    plt.legend(title='Amplification Type', loc='upper right')
    output_plot = os.path.join(output_folder, f'MVD_{filename_prefix}_all_key_pairs.png')
    plt.savefig(output_plot)
    plt.close()

    print('Boxplot for all key pairs was saved')

###############################################################################
## Allele dropout
def calculate_ado(col):
    return (col == 0).mean() * 100

## plot ADO
def plot_ado(csv_file, key1, key2, output_folder):
    # Load data
    dp_vaf_data = pd.read_csv(csv_file, sep='\t')
    filename_prefix = os.path.splitext(os.path.basename(csv_file))[0]
    output_plot = os.path.join(output_folder, f'ADO_{filename_prefix}.png')
    
    # Filter MDA and PTA VAF columns
    mda_vaf = dp_vaf_data.filter(like=f'VAF_{key1}')
    pta_vaf = dp_vaf_data.filter(like=f'VAF_{key2}')
    
    # Calculate ADO percentages
    mda_ado = mda_vaf.apply(calculate_ado)
    pta_ado = pta_vaf.apply(calculate_ado)
    
    # Prepare data for plotting
    ado_data = pd.concat([
        pd.DataFrame({'Sample': mda_ado.index, 'ado': mda_ado.values, 'Sample Type': 'MDA'}),
        pd.DataFrame({'Sample': pta_ado.index, 'ado': pta_ado.values, 'Sample Type': 'PTA'})
    ])

    # Plotting
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Sample Type', y='ado', data=ado_data, palette='husl')
    plt.xlabel('Simulations')
    plt.ylabel('ADO (%)')
    plt.title('Box Plot of ADO for MDA and PTA Samples')
    plt.legend(title='Sample Type', loc='upper right')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.show()
    print('Boxplot was saved to:', output_plot)

###############################################################################
## plotting the depth
def depth_plot(csv_file, output_folder):
    # Load data
    dp_vaf_data = pd.read_csv(csv_file, sep='\t')
    filename_prefix = os.path.splitext(os.path.basename(csv_file))[0]
    output_plot = os.path.join(output_folder, f'DP_{filename_prefix}.png')
    
    # Filter and reshape depth data for MDA and PTA
    depth_data = pd.concat([
        dp_vaf_data.filter(like='DP_MDA').melt(var_name='Sample', value_name='Depth').assign(Sample_Type='MDA'),
        dp_vaf_data.filter(like='DP_PTA').melt(var_name='Sample', value_name='Depth').assign(Sample_Type='PTA')
    ])
    
    # Create box plot
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Sample', y='Depth', hue='Sample_Type', data=depth_data, palette='husl')
    plt.xlabel('Sample')
    plt.ylabel('Depth')
    plt.title('Box Plot of Depth for MDA and PTA Samples')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_plot)
    plt.show()  # Optional: show the plot after saving
    print('Boxplot was saved to:', output_plot)
    
###############################################################################
