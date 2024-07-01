import sys
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt

def calculate_statistic(row):
    try:
        vaf_mda_values = row.filter(like='MDA').filter(like='VAF').values 
        vaf_pta_values = row.filter(like='PTA').filter(like='VAF').values
        vaf_mda_values = np.abs(vaf_mda_values - 0.5)
        vaf_pta_values = np.abs(vaf_pta_values - 0.5)
        mda_statistic = np.median(vaf_mda_values)
        pta_statistic = np.median(vaf_pta_values)
        return pd.Series({'MDA_Statistic': mda_statistic, 'PTA_Statistic': pta_statistic})
    except Exception as e:
        print("Error:", e)
        return pd.Series({'MDA_Statistic': np.nan, 'PTA_Statistic': np.nan})
    
    
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input_file.csv")
        sys.exit(1)
    input_file = sys.argv[1]
    merged_data = pd.read_csv(input_file, sep= '\t')
    statistics_df = merged_data.apply(calculate_statistic, axis=1)
    statistics_df.dropna(inplace=True)

    statistic, p_value = mannwhitneyu(statistics_df['MDA_Statistic'], statistics_df['PTA_Statistic'])

    print("Mann-Whitney U test:")
    print(f"Statistic: {statistic}")
    print(f"P-value: {p_value}")

    alpha = 0.05
    if p_value < alpha:
        print("Reject null hypothesis: There is a significant difference between the groups.")
    else:
        print("Fail to reject null hypothesis: There is no significant difference between the groups.")

    results_df = pd.DataFrame({'Statistic': [statistic], 'P-value': [p_value]})

    # Add median columns to the results DataFrame
    
    results_df['MDA_Median'] = statistics_df['MDA_Statistic']
    results_df['PTA_Median'] = statistics_df['PTA_Statistic']

    # Save the results DataFrame to CSV with median columns
    results_df.to_csv("mann_whitney_results.csv", index=False)
    print("Results saved to 'mann_whitney_results.csv'.")

    # Extract the medians for MDA and PTA groups
    group_medians = statistics_df[['MDA_Statistic', 'PTA_Statistic']]

    # Plot the distribution of medians using Seaborn
    sns.boxplot(data=group_medians)
    plt.xlabel('Group')
    plt.ylabel('Median')
    plt.title('Distribution of Medians for MDA and PTA Groups')

    # Save the plot to a file
    plt.savefig("medians_distribution_plot.png")

    # Extract the MDA and PTA depth values
    mda_depth = merged_data.filter(like='DP_MDA')
    pta_depth = merged_data.filter(like='DP_PTA')

    # Reshape the data to long format
    mda_depth = mda_depth.melt(var_name='Sample', value_name='Depth')
    mda_depth['Sample Type'] = 'MDA'
    pta_depth = pta_depth.melt(var_name='Sample', value_name='Depth')
    pta_depth['Sample Type'] = 'PTA'

    # Concatenate MDA and PTA data
    depth_data = pd.concat([mda_depth, pta_depth])

    # Plot the box plot
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Sample', y='Depth', hue='Sample Type', data=depth_data, palette='Set2')
    plt.xlabel('Sample')
    plt.ylabel('Depth')
    plt.title('Box Plot of Depth for MDA and PTA Samples')
    plt.legend(title='Sample Type')
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save the plot to a file
    plt.savefig("depth_boxplot.png")
    
    statistics_df['Difference'] = statistics_df['MDA_Statistic'] - statistics_df['PTA_Statistic']

    # Plot the boxplot of differences
    plt.figure(figsize=(8, 6))
    sns.boxplot(y='Difference', data=statistics_df)
    plt.xlabel('Difference')
    plt.ylabel('Values')
    plt.title('Boxplot of Differences between MDA and PTA Statistics')
    plt.savefig("differences_boxplot.png")
    plt.show()