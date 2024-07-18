import argparse
import multiprocessing
import time
import sys
import os
import logging
import numpy as np
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)
from SinglCellSim.AmpSim.MDAstate import MDASimulation, subsetAmpliconSaveToFASTA
from SinglCellSim.configs.configFunctions import read_config, read_fasta

logging.basicConfig(level=logging.DEBUG)  # Set DEBUG level to capture all messages

def run_simulation(sim_index, args, config_params):
    try:
        logging.debug(f"Running simulation {sim_index} with PID {os.getpid()}")
        logging.debug(f"Config params at start of simulation {sim_index}: {config_params}")
        
        np.random.seed((os.getpid() * int(time.time())) % 123456789)
        output_folder = os.path.join(args.output_base, f"sim{sim_index}")
        os.makedirs(output_folder, exist_ok=True)
        
        Theta = int(config_params['theta'])
        Gamma = int(config_params['gamma'])
        DNACoef = int(config_params['dnacoef'])
        lMin = int(config_params['lmin'])
        lMax = int(config_params['lmax'])
        Lambda = float(config_params['lambda'])
        delta_t = float(config_params['delta_t'])
        beta = float(config_params['beta'])
        exclude = config_params['exclude']

        start_time = time.time()

        # Read patSeq and matSeq from FASTA files
        patSeq_data = read_fasta(args.patSeq_file)
        matSeq_data = read_fasta(args.matSeq_file)
        template = args.tem
        
        amplicons = MDASimulation(patSeq_data, matSeq_data, Theta=Theta, Gamma=Gamma, DNACoef=DNACoef,
                                  lMin=lMin, lMax=lMax, Lambda=Lambda, delta_t=delta_t, beta=beta, 
                                  exclude=exclude, template=template, output_folder=output_folder)
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        logging.debug(f"Elapsed time for simulation {sim_index}: {elapsed_time} seconds")

        subsetAmpliconSaveToFASTA(amplicons, patSeq_data, matSeq_data, output_folder=output_folder)
    except Exception as e:
        logging.error(f"Error in simulation {sim_index}: {e}")

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

    logging.debug(f"Configuration parameters: {config_params}")

    try:
        # Temporarily run simulations sequentially for debugging
        for sim_index in range(1, args.num_simulations + 1):
            run_simulation(sim_index, args, config_params)

        # Uncomment below for multiprocessing after debugging
        # with multiprocessing.Pool(processes=args.num_cores) as pool:
        #     pool.starmap(run_simulation, [(sim_index, args, config_params) for sim_index in range(1, args.num_simulations + 1)])
    except Exception as e:
        logging.error(f"Error in multiprocessing: {e}")

    logging.debug("All simulations completed.")

if __name__ == "__main__":
    main()
