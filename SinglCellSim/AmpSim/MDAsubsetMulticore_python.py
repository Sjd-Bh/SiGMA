import argparse
import multiprocessing
import time
import sys
import os
import numpy as np
np.random.seed(42)  # Set a specific seed value, such as 42
import random
random.seed(42)
sys.path.append('/home/bahonar/simulation/SingleCellSim')
from AmpSim.MDAstate import MDASimulation, subsetAmpliconSaveToFASTA
from configs.configFunctions import read_config, read_fasta

def run_simulation(sim_index, args, config_params):
    np.random.seed((os.getpid() * int(time.time())) % 123456789)
    output_folder = os.path.join(args.output_base, f"sim{sim_index}")
    os.makedirs(output_folder, exist_ok=True)
    Theta = config_params['theta']
    Gamma = config_params['gamma']
    DNACoef = config_params['dnacoef']
    lMin = config_params['lmin']
    lMax = config_params['lmax']
    Lambda = config_params['lambda']
    delta_t = config_params['delta_t']
    beta = config_params['beta']
    exclude = config_params['exclude']

    start_time = time.time()

    # Read patSeq and matSeq from FASTA files
    patSeq_data = read_fasta(args.patSeq_file)
    matSeq_data = read_fasta(args.matSeq_file)
    template = args.tem
    
    amplicons = MDASimulation(patSeq_data, matSeq_data, Theta=Theta, Gamma=Gamma, DNACoef=DNACoef,
                              lMin=lMin, lMax=lMax, Lambda=Lambda, delta_t=delta_t, beta=beta, 
                              exclude=exclude, template = template,output_folder=output_folder)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for simulation {sim_index}: {elapsed_time} seconds")

    subsetAmpliconSaveToFASTA(amplicons, patSeq_data, matSeq_data, output_folder=output_folder)

def main():
    parser = argparse.ArgumentParser(description="Run MDA simulations and save subsets to FASTA")
    parser.add_argument("--num_simulations", type=int, default=1, help="Number of simulations to run")
    parser.add_argument("--config_file", type=str, default="MDAsim.ini", help="Path to config.ini file")
    parser.add_argument("--output_base", type=str, default="output", help="Base path to the output folder")
    parser.add_argument("--patSeq_file", type=str, help="Path to paternal genome FASTA file")
    parser.add_argument("--matSeq_file", type=str, help="Path to maternal genome FASTA file")
    parser.add_argument("--num_cores", type=int, default=multiprocessing.cpu_count(), help="Number of CPU cores to use")
    parser.add_argument("--log_file", type=str, default="output.log", help="Path to the log file")
    parser.add_argument("--resume", action="store_true", help="Resume from the last checkpoint")
    parser.add_argument("--depth", type=int, help="the desired depth")
    parser.add_argument("--tem", action="store_true", help="amplification template from template or not")
    args = parser.parse_args()

    # Read parameters from the configuration file
    config_params = read_config(args.config_file, "Simulation")

    # Convert string values to appropriate types
    config_params = {key: int(value) if value.isdigit() else float(value) if '.' in value else value
                     for key, value in config_params.items()}
    print("Keys in config_params:", config_params.keys())
    np.random.seed((os.getpid() * int(time.time())) % 123456789)
    with multiprocessing.Pool(processes=args.num_cores) as pool:
        pool.starmap(run_simulation, [(sim_index, args, config_params) for sim_index in range(1, args.num_simulations + 1)])

    print("All simulations completed.")

if __name__ == "__main__":
    main()