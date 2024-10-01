import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# def corFunc(vafCol, posDiff, min_val, max_val):
#     indexes = np.where((posDiff > min_val) & (posDiff < max_val))[0]
#     return np.corrcoef(vafCol[indexes], vafCol[indexes - 1])[0, 1]

# def pos_vaf_separator(DPvafData,vafData ,colname_POS,colname_VAF):
#     pos = DPvafData[colname_POS].values
#     vaf = vafData[colname_VAF].values
#     return pos, vaf  

##############################################################
def calculate_correlations(posDiff,data, min_val, max_val):
    correlations = []
    for col in data.columns:
        selected_col = data[col]
        indexes = np.where((posDiff > min_val) & (posDiff < max_val))[0]
        correlation = np.corrcoef(selected_col[indexes], selected_col[indexes - 1])[0, 1]
        correlations.append(correlation)
    return correlations
##############################################################
def MDA_PTA_VAF_DP_separation(DPvafData,key1, key2):
    # Filter columns based on conditions and extract their names
    VAF_MDA_cols = DPvafData.filter(like=key1).filter(like='VAF').columns
    VAF_PTA_cols = DPvafData.filter(like=key2).filter(like='VAF').columns 
    DP_MDA_cols = DPvafData.filter(like=key1).filter(like='DP').columns 
    DP_PTA_cols = DPvafData.filter(like=key2).filter(like='DP').columns
    
    # Extract data based on filtered column names
    VAF_MDA = DPvafData[VAF_MDA_cols]
    VAF_PTA = DPvafData[VAF_PTA_cols]
    DP_MDA = DPvafData[DP_MDA_cols]
    DP_PTA = DPvafData[DP_PTA_cols]
    posCol_diff = DPvafData['POS'].diff()
    posCol_diff = posCol_diff.fillna(0)
    
    # Return the four sets of data as a tuple
    return VAF_MDA, VAF_PTA, DP_MDA, DP_PTA, posCol_diff
############################################################## 
def plot_boxplot(data, labels, min_values, max_values, title, filename):
    plt.figure(figsize=(10, 5))
    
    if len(data) == 1:
        plt.boxplot(data[0], positions=[1], labels=labels, patch_artist=True)
        plt.gca().artists[0].set_facecolor('blue')  # Change 'blue' to desired color
    else:
        positions = np.arange(1, len(data) * 2 + 1, 2)
        plt.boxplot(data, positions=positions, labels=labels, patch_artist=True)
        
        # Adding color to boxes
        num_boxes = len(data)
        colors = ['#FF00FF', '#9A0EEA'] * (num_boxes // 2)
        for box, color in zip(plt.gca().artists, colors):
            box.set_facecolor(color)
    
    # Set x-axis labels with min and max values
    for i, (min_val, max_val) in enumerate(zip(min_values, max_values)):
        label = f'{labels[i]}\nMin: {min_val}\nMax: {max_val}'
        plt.text(2 * i + 1, -0.1, label, ha='center')
    
    plt.title(title)
    plt.xlabel('DataFrames')
    plt.ylabel('Correlation Coefficient')
    plt.grid(axis='y')
    # plt.tight_layout()
    # plt.savefig(filename)
    # print(f'file saved in {filename}')
    plt.close()
#############################################################
def plot_boxplots_with_ranges(vafMDA, vafPTA, dpMDA, dpPTA, min_values, max_values, filename_prefix, output_folder,posDiff):
    # Define lists to store correlations for VAF and DP
    correlations_vaf = {'MDA': [], 'PTA': []}
    correlations_dp = {'MDA': [], 'PTA': []}
    
    # Iterate over min and max value pairs
    for min_val, max_val in zip(min_values, max_values):
        # Calculate correlations for VAF and DP based on current min and max values
        correlations_vaf['MDA'].append(calculate_correlations(posDiff,vafMDA, min_val, max_val))
        correlations_vaf['PTA'].append(calculate_correlations(posDiff,vafPTA, min_val, max_val))
        correlations_dp['MDA'].append(calculate_correlations(posDiff,dpMDA, min_val, max_val))
        correlations_dp['PTA'].append(calculate_correlations(posDiff,dpPTA, min_val, max_val))
        
    # Create labels for box plots (alternating between MDA and PTA)
    num_pairs = len(min_values)
    labels_vaf = ['VAF MDA'] * num_pairs + ['VAF PTA'] * num_pairs
    labels_dp = ['DP MDA'] * num_pairs + ['DP PTA'] * num_pairs

    # Flatten the correlations lists for plotting
    correlations_vaf_flat = correlations_vaf['MDA'] + correlations_vaf['PTA']
    correlations_dp_flat = correlations_dp['MDA'] + correlations_dp['PTA']

    outputPlotVAF = os.path.join(output_folder, f'VAF_cor_{filename_prefix}.png')
    outputPlotDP = os.path.join(output_folder, f'DO_cor_{filename_prefix}.png')
    # Plot box plots for VAF and DP with the prefix based on the CSV file name
    plot_boxplot(correlations_vaf_flat, labels_vaf,min_values, max_values, 'Correlations for VAF', outputPlotVAF)
    plot_boxplot(correlations_dp_flat, labels_dp,min_values, max_values, 'Correlations for DP', outputPlotDP)

# Example usage:
# plot_boxplots_with_ranges(vafMDA, vafPTA, dpMDA, dpPTA, min_values, max_values, 'DPvafData.csv')
#############################################################
def ADO(col):
    zeros_count = (col == 0).sum()
    total_count = len(col)
    percentage_zeros = (zeros_count / total_count) * 100
    return percentage_zeros

#############################################################

# def calculate_statistic(row):
#     try:
#         vaf_mda_values = row.filter(like='MDA').filter(like='VAF').values 
#         vaf_pta_values = row.filter(like='PTA').filter(like='VAF').values
#         vaf_mda_values = np.abs(vaf_mda_values - 0.5)
#         vaf_pta_values = np.abs(vaf_pta_values - 0.5)
#         mda_statistic = np.median(vaf_mda_values)
#         pta_statistic = np.median(vaf_pta_values)
#         return pd.Series({'MDA_Statistic': mda_statistic, 'PTA_Statistic': pta_statistic})
#     except Exception as e:
#         print("Error:", e)
#         return pd.Series({'MDA_Statistic': np.nan, 'PTA_Statistic': np.nan})
    
    
    
def calculate_statistic(row,key1,key2):
    vaf_Allow_values = row.filter(like=f'VAF_{key1}').values
    vaf_notAllow_values = row.filter(like=f'VAF_{key2}').values 
    vaf_Allow_values = np.abs(vaf_Allow_values - 0.5)
    vaf_notAllow_values = np.abs(vaf_notAllow_values - 0.5)
    Allow_statistic = np.median(vaf_Allow_values)
    notAllow_statistic = np.median(vaf_notAllow_values)
    return pd.Series({'Statistic1': Allow_statistic, 'Statistic2': notAllow_statistic})
#############################################################

def wilcoxonMDAPTA(csvFiles, key1,key2,output_folder):
    mergedData = pd.DataFrame()
    for csvFile in csvFiles:
        DPvafData = pd.read_csv(csvFile, sep= '\t')
        filename_prefix = os.path.splitext(os.path.basename(csvFile))[0]
        # keys = [(key1, key2), (key1, key3), (key1,key4), (key3, key4)]
        keys = [ (key1, key2)]
        for keyPair in keys:
            
            statistics_df = DPvafData.apply(calculate_statistic,args = keyPair,  axis=1)
            statistics_df.dropna(inplace=True)
            statistics_df['Difference'] = statistics_df['Statistic1'] - statistics_df['Statistic2']
            statType = f'{keyPair}'
            statistics_df['statType'] = statType
            statistics_df['amp'] = filename_prefix
            # mergedData['amp'] = pd.Series([filename_prefix] * len(mergedData))
            Data = statistics_df[['Difference','statType','amp']]
            mergedData = pd.concat([mergedData,Data],ignore_index = True)
        
    outputCSV = os.path.join(output_folder, f'wilcoxon_{filename_prefix}.csv')
    outputPlot = os.path.join(output_folder, f'wilcoxon_{filename_prefix}.png')
    mergedData.to_csv(outputCSV, index=False, sep='\t')
    
    # Plot the boxplot of differences
    plt.figure(figsize=(8, 6))
    sns.set_theme(style="ticks", palette="pastel")
    # ax = sns.boxplot(x='amp', y='Difference', hue='statType', data=mergedData, palette="Set1", width=0.6, whis=1.5, showfliers=False)
    ax = sns.boxplot(x='amp', y='Difference', data=mergedData, palette="Set1", width=0.6, whis=1.5, showfliers=False)
    
    # Calculate and annotate mean for each boxplot
    # group_means = mergedData.groupby(['amp', 'statType'])['Difference'].mean()
    # for i, (label, mean_value) in enumerate(group_means.items()):
    #     ax.text(i, mean_value, f'Mean: {mean_value:.2f}', fontsize=10, color='black', ha='center', va='bottom')
    plt.xlabel('wilcoxn')
    plt.xticks(rotation=45)
    plt.ylabel('Difference')
    plt.title('')
    plt.savefig(outputPlot)
    plt.show()
    print('boxplot was saved')
##################################################################
def depthPlot(csvFile, key1,key2,output_folder ):
    DPvafData = pd.read_csv(csvFile, sep= '\t')
    filename_prefix = os.path.splitext(os.path.basename(csvFile))[0]
    outputPlot = os.path.join(output_folder, f'DP_{filename_prefix}.png')
    
    #### depth
    mda_depth = DPvafData.filter(like= 'DP_MDA')
    pta_depth = DPvafData.filter(like='DP_PTA')

    # Reshape the data to long format
    mda_depth = mda_depth.melt(var_name='Sample', value_name='Depth')
    mda_depth['Sample Type'] = 'MDA'
    pta_depth = pta_depth.melt(var_name='Sample', value_name='Depth')
    pta_depth['Sample Type'] = 'PTA'

    # Concatenate MDA and PTA data
    depth_data = pd.concat([mda_depth, pta_depth])
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Sample', y='Depth', hue='Sample Type', data=depth_data, palette='husl')
    plt.xlabel('Sample')
    plt.ylabel('Depth')
    plt.title('Box Plot of Depth for MDA and PTA Samples')
    plt.legend(title='Sample Type')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(outputPlot)

###################################################################
def corPlot(patSNPfile, csvFile,min_values, max_values,key1, key4, output_folder):
    DPvafData = pd.read_csv(csvFile, sep= '\t')
    filename_prefix = os.path.splitext(os.path.basename(csvFile))[0]
    #phasing
    patSNP = pd.read_csv(patSNPfile, sep='\t', comment='#', header=None)
    column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    patSNP.columns = column_names
    snp_positions = sorted(patSNP['POS'])
    for pos in snp_positions:
        if pos in DPvafData['POS'].values:
            for col in DPvafData.columns:
                if 'VAF' in col:
                    DPvafData.loc[DPvafData['POS'] == pos, col] = 1 - DPvafData.loc[DPvafData['POS'] == pos, col]
        
        
    # boxplot for correlations
    vafMDA, vafPTA, dpMDA, dpPTA, posDiff = MDA_PTA_VAF_DP_separation(DPvafData, key1, key4)
    corvaf = pd.DataFrame()
    cordp = pd.DataFrame()
    
    for min_val, max_val in zip(min_values, max_values):
        cor_vaf_MDA = calculate_correlations(posDiff,vafMDA, min_val, max_val)
        cor_vaf_PTA = calculate_correlations(posDiff,vafPTA, min_val, max_val)
        cor_dp_MDA = calculate_correlations(posDiff,dpMDA, min_val, max_val)
        cor_dp_PTA = calculate_correlations(posDiff,dpPTA, min_val, max_val)
        
        corvafMDA = pd.DataFrame(cor_vaf_MDA,columns=['cor'])
        corvafPTA = pd.DataFrame(cor_vaf_PTA,columns=['cor'])
        cordpMDA = pd.DataFrame(cor_dp_MDA,columns=['cor'])
        cordpPTA = pd.DataFrame(cor_dp_PTA,columns=['cor'])
        corvafMDA['range'] = f'{min_val}_{max_val}'
        corvafMDA['type'] = 'MDA'
        corvafPTA['range'] = f'{min_val}_{max_val}'
        corvafPTA['type'] = 'PTA'
        cordpMDA['range'] = f'{min_val}_{max_val}'
        cordpMDA['type'] = 'MDA'
        cordpPTA['range'] = f'{min_val}_{max_val}'
        cordpPTA['type'] = 'PTA'
        
        corvaf = pd.concat([corvaf,corvafMDA],ignore_index = True)
        corvaf = pd.concat([corvaf,corvafPTA],ignore_index = True)
        
        cordp = pd.concat([cordp,cordpMDA],ignore_index = True)
        cordp = pd.concat([cordp,cordpPTA],ignore_index = True)
        
    # plt.figure(figsize=(8, 6))
    # sns.boxplot(x='type', y='cor', hue='range', data=corvaf, palette="Set1", width=0.6, whis=1.5, showfliers=False)
    # # ax = sns.boxplot(x='amp', y='Difference', data=mergedData, palette="Set1", width=0.6, whis=1.5, showfliers=False)
    
    # # Calculate and annotate mean for each boxplot
    # # group_means = mergedData.groupby(['amp', 'statType'])['Difference'].mean()
    # # for i, (label, mean_value) in enumerate(group_means.items()):
    # #     ax.text(i, mean_value, f'Mean: {mean_value:.2f}', fontsize=10, color='black', ha='center', va='bottom')
    # plt.xlabel('correlations across repeat simulations')
    # plt.xticks(rotation=45)
    # plt.ylabel('correlations')
    # plt.title('correlations of vaf based on physical distances')
    # # plt.savefig(outputPlot)
    # plt.show()
    # print('boxplot was saved')   
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='type', y='cor', hue='range', data=corvaf, palette="Set1", width=0.8, whis=1.5, showfliers=False)
    # Increase the width parameter from 0.6 to 0.8 or adjust as needed
    # You can play around with different values to achieve the desired spacing
    plt.xlabel('correlations across repeat simulations')
    plt.xticks(rotation=45)
    plt.ylabel('correlations')
    plt.title('correlations of vaf based on physical distances')
    plt.show()
    print('boxplot was saved')
        
     
        
        
    
###################################################################
def adoPlot(csvFile,key1,key2,output_folder):
    DPvafData = pd.read_csv(csvFile, sep= '\t')
    filename_prefix = os.path.splitext(os.path.basename(csvFile))[0]
    outputPlot = os.path.join(output_folder, f'ADO_{filename_prefix}.png')
    mda_vaf = DPvafData.filter(like= f'VAF_{key1}')
    pta_vaf = DPvafData.filter(like= f'VAF_{key2}')
    mdaADO = pd.DataFrame(mda_vaf.apply(ADO,  axis=0))
    ptaADO = pd.DataFrame(pta_vaf.apply(ADO,  axis=0))
    mdaADO = mdaADO.melt(var_name='Sample', value_name='ado')
    mdaADO['Sample Type'] = 'MDA'
    ptaADO = ptaADO.melt(var_name='Sample', value_name='ado')
    ptaADO['Sample Type'] = 'PTA'
    
    ADOdata = pd.concat([mdaADO, ptaADO])
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Sample Type', y='ado', data=ADOdata, palette='husl')
    plt.xlabel('simulations')
    plt.ylabel('ADO')
    plt.title('Box Plot of ADO for MDA and PTA Samples')
    plt.legend(title='Sample Type')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(outputPlot)
    
#####################################################################
    
  