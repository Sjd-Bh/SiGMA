import argparse
#from processDP_VAF_CSVtoPlot import depthPlot, corPlot, adoPlot, MVD_Diff_MP
import sys
import os

# Get the absolute path to the project directory based on the current script's location
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
print(project_path)
sys.path.append(project_path)

from SinglCellSim.configs.configFunctions import read_config
from evaluation.comparePlots.simSCVtoPlot import MDA_PTA_VAF_DP_separation,calculate_ado,calculate_correlations,calculate_mvd, cor_plot,depth_plot,mvd_diff_mp,plot_ado  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SNP calling on single cells')
    parser.add_argument('--csv_files', nargs='+', type=str, help='List of paths to single-cell data files')
    parser.add_argument('--output_folder', type=str, help='output file path and name')
    parser.add_argument('--key1', type=str, help='MDA or COCAllow')
    parser.add_argument('--key2', type=str, help='PTA or COCnotAllow')
    parser.add_argument("--config_file", type=str, default="mvdConfig.ini", help="Path to config.ini file")
    args = parser.parse_args()
    csvFiles = args.csv_files
    key1 = args.key1
    key2 = args.key2
    output_folder = args.output_folder
    config_params = read_config(args.config_file, "Plot")
    min_values = config_params['min_values']
    max_values = config_params['max_values']

    for csvFile in csvFiles:
        mvd_diff_mp(csvFile, key1, key2, output_folder)
        depth_plot(csvFile, key1, key2, output_folder)
        cor_plot(csvFile, min_values, max_values, key1, key2, output_folder)
        plot_ado(csvFile, key1, key2, output_folder)
        
    print("All plots of CSV files are saved")
