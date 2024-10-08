import argparse
#from processDP_VAF_CSVtoPlot import depthPlot, corPlot, adoPlot, MVD_Diff_MP
import sys
import os

# Get the absolute path to the project directory based on the current script's location
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

from evaluation.configs.configFunctions import read_config
from evaluation.comparePlots.simSCVtoPlot import MDA_PTA_VAF_DP_separation,calculate_ado,calculate_correlations,calculate_mvd, cor_plot,depth_plot,mvd_diff_mp,plot_ado  

import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SNP calling on single cells')
    parser.add_argument('--csv_files', nargs='+', type=str, help='List of paths to single-cell data files or directories containing them')
    parser.add_argument('--output_folder', type=str, help='output file path and name')
    parser.add_argument('--key1', type=str, help='MDA or COCAllow')
    parser.add_argument('--key2', type=str, help='PTA or COCnotAllow')
    parser.add_argument("--config_file", type=str, default="mvdConfig.ini", help="Path to config.ini file")
    args = parser.parse_args()
    
    # Collect all CSV files from provided directories or file paths
    csv_files = []
    for path in args.csv_files:
        if os.path.isdir(path):
            # Recursively find all CSV files in the directory
            csv_files.extend(glob.glob(os.path.join(path, "*.csv")))
        else:
            # Add file directly if it's a file
            csv_files.append(path)

    key1 = args.key1
    key2 = args.key2
    output_folder = args.output_folder
    config_params = read_config(args.config_file, "Plot")
    min_values = config_params['min_values']
    max_values = config_params['max_values']

    for csvFile in csv_files:
        mvd_diff_mp(csvFile, key1, key2, output_folder)
        depth_plot(csvFile, key1, key2, output_folder)
        cor_plot(csvFile, min_values, max_values, key1, key2, output_folder)
        plot_ado(csvFile, key1, key2, output_folder)
        
    print("All plots of CSV files are saved")