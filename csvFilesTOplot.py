import argparse
from processDP_VAF_CSVtoPlot import depthPlot, corPlot, adoPlot, wilcoxonMDAPTA
import sys
sys.path.append('/home/bahonar/simulation/SingleCellSim')
from configs.configFunctions import read_config

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SNP calling on single cells')
    parser.add_argument('--csv_files', nargs='+', type=str, help='List of paths to single-cell data files')
    parser.add_argument('--output_folder', type=str, help='output file path and name')
    parser.add_argument('--key1', type=str, help='MDA or COCAllow')
    parser.add_argument('--key2', type=str, help='PTA or COCnotAllow')
    parser.add_argument("--config_file", type=str, default="wilcoxonConfig.ini", help="Path to config.ini file")
    args = parser.parse_args()
    csvFiles = args.csv_files
    key1 = args.key1
    key2 = args.key2
    output_folder = args.output_folder
    config_params = read_config(args.config_file, "Plot")
    min_values = config_params['min_values']
    max_values = config_params['max_values']
 
    for csvFile in csvFiles:
        wilcoxonMDAPTA(csvFile, key1,key2, output_folder)
        depthPlot(csvFile, key1,key2,output_folder )
        corPlot(csvFile,min_values, max_values,key1, key2, output_folder)
        adoPlot(csvFile,key1,key2,output_folder)
        
    print("all plots of csv files are saved")    