import argparse
#from processDP_VAF_CSVtoPlot import depthPlot, corPlot, adoPlot, MVD_Diff_MP
import sys
import os
import ast

# Get the absolute path to the project directory based on the current script's location
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

from evaluation.configs.configFunctions import read_config
from evaluation.realDataAnalysis.csvToPlot import MDA_PTA_VAF_DP_separation,calculate_ado,calculate_correlations,calculate_mvd, cor_plot,depth_plot,mvd_diff_mp,plot_ado  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SNP calling on single cells')
    parser.add_argument('--csv_files', nargs='+', type=str, help='List of paths to single-cell data files')
    parser.add_argument('--output_folder', type=str, help='output file path and name')
    parser.add_argument('--key1', type=str, help='MDA or COCAllow')
    parser.add_argument('--key2', type=str, help='PTA or COCnotAllow')
    parser.add_argument("--config_file", type=str, default="mvdConfig.ini", help="Path to config.ini file")
    # parser.add_argument("--pat-SNP", type=str, default="", help="paternal SNPs for phasing")
    
    args = parser.parse_args()
    csvFiles = args.csv_files
    key1 = args.key1
    key2 = args.key2
    output_folder = args.output_folder
    config_params = read_config(args.config_file, "Plot")
    # pat_SNP_file = args.pat_SNP
    min_values_str = config_params['min_values']
    max_values_str = config_params['max_values']
    
    # Convert the string lists to Python lists of integers using ast.literal_eval
    min_values = ast.literal_eval(min_values_str)
    max_values = ast.literal_eval(max_values_str)
    # Ensure they are numeric (just in case)
    min_values = [int(value) for value in min_values]
    max_values = [int(value) for value in max_values]
    
    mvd_diff_mp(csvFiles, key1, key2, output_folder)
    
    for csvFile in csvFiles:
        depth_plot(csvFile, output_folder)
        cor_plot(csvFile, min_values, max_values, key1, key2, output_folder)
        plot_ado(csvFile, key1, key2, output_folder)
        
    print("All plots of CSV files are saved")
