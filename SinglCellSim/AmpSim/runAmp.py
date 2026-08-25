import argparse
import os
import sys
import time
import multiprocessing
import random
import traceback
import pickle
import numpy as np

# --------------------------------------------------
# Project setup
# --------------------------------------------------
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(
    os.path.join(current_dir, os.pardir, os.pardir)
)
sys.path.insert(0, project_root)

from SinglCellSim.configs.configFunctions import read_config, read_fasta
from SinglCellSim.AmpSim.AmpFunction import (
    MDASimulation,
    subsetAmpliconSaveToFASTA,
    load_cnvs_from_bed
)

# --------------------------------------------------
# Single simulation runner
# --------------------------------------------------
def run_simulation(sim_id, args, config):
    """
    Run one CNV-aware two-phase MDA simulation:
    Phase 1: initiation on parental strands
    Phase 2: exponential amplification on amplicons
    """
    try:
        # ------------------------------
        # Reproducible seed
        # ------------------------------
        seed = (int(time.time() * 1e6) % (2**32)) + sim_id
        np.random.seed(seed)
        random.seed(seed)

        print(f"[Sim {sim_id}] Seed = {seed}")

        # ------------------------------
        # Load reference genomes
        # ------------------------------
        pat_seq = str(read_fasta(args.patSeq_file))
        mat_seq = str(read_fasta(args.matSeq_file))

        # ------------------------------
        # Simulation parameters
        # ------------------------------
        lMin        = int(config["lmin"])
        lMax        = int(config["lmax"])
        delta_t     = float(config["delta_t"])
        beta        = float(config["beta"])
        total_time  = float(config["total_time"])
        lambda_init = float(config["lambda_init"])
        lambda_exp  = float(config["lambda_exp"])
        Theta_init  = int(config["theta_init"])
        Theta_exp   = int(config["theta_exp"])
        t_switch    = int(config["t_switch"])
        
        # ------------------------------
        # Load CNVs
        # ------------------------------
        cnvs, cnvs_for_mda = load_cnvs_from_bed(args.cnv_bed_file)
        
        # ------------------------------
        # Output directory
        # ------------------------------
        output_folder = os.path.join(args.output_base, f"sim{sim_id}")
        os.makedirs(output_folder, exist_ok=True)
        
        # ------------------------------
        # Run two-phase MDA simulation
        # ------------------------------
        sim_start = time.time()

        (
            amplicons,        # final amplicon list
            cycle_stats,      # per-cycle statistics
            P_count,          # parental P usage
            M_count,          # parental M usage
            cnvs_used,        # CNVs actually applied
            bp_P,             # Paternal breakpoints (replaces segments_P)
            bp_M              # Maternal breakpoints (replaces segments_M)
        ) = MDASimulation(
            patSeq=pat_seq,
            matSeq=mat_seq,
            lMin=lMin,
            lMax=lMax,
            delta_t=delta_t,
            beta=beta,
            CNVs=cnvs_for_mda,
            total_time=total_time,
            t_switch=t_switch,
            lambda_init=lambda_init,
            lambda_exp=lambda_exp,
            Theta_init=Theta_init,
            Theta_exp=Theta_exp,
            output_folder=output_folder
        )

        sim_end = time.time()
        total_runtime = sim_end - sim_start

        # ------------------------------
        # Runtime statistics
        # ------------------------------
        n_cycles = len(cycle_stats)
        avg_cycle_time = total_runtime / n_cycles if n_cycles > 0 else 0.0

        with open(os.path.join(output_folder, "cycle_runtime.log"), "w") as f:
            f.write("#cycle\tavg_cycle_sec\tcumulative_sec\n")
            cumulative = 0.0
            for cycle, phase, _, _ in cycle_stats:
                cumulative += avg_cycle_time
                f.write(f"{cycle}\t{avg_cycle_time:.6f}\t{cumulative:.6f}\n")

        # ------------------------------
        # Save FASTA output
        # ------------------------------
        # The new approach writes directly from the reference string
        # since errors and bounds are native to reference coordinates.
        subsetAmpliconSaveToFASTA(
            amplicons=amplicons,
            refSeq_P=pat_seq,
            refSeq_M=mat_seq,
            lMin=lMin,
            output_folder=output_folder
        )        
        
        # ------------------------------
        # Save cycle statistics
        # ------------------------------
        with open(os.path.join(output_folder, "amplicon_stats.tsv"), "w") as f:
            f.write("cycle\tnew_amplicons\ttotal_length_bp\n")
            for cycle, phase, n_new, total_len in cycle_stats:
                f.write(f"{cycle}\t{n_new}\t{total_len}\n")

        # ------------------------------
        # Save full simulation state
        # ------------------------------
        with open(os.path.join(output_folder, "amplicons.pkl"), "wb") as f:
            pickle.dump(amplicons, f)

        print(
            f"[Sim {sim_id}] Finished in {total_runtime:.2f}s "
            f"({n_cycles} cycles, ~{avg_cycle_time:.4f}s/cycle)"
        )

    except Exception as e:
        print(f"[Sim {sim_id}] FAILED: {e}")
        traceback.print_exc()

# --------------------------------------------------
# Multiprocessing wrapper
# --------------------------------------------------
def run_all_simulations(args, config):
    print(
        f"\nRunning {args.num_simulations} simulations "
        f"on {args.num_cores} cores\n"
    )

    with multiprocessing.Pool(processes=args.num_cores) as pool:
        pool.starmap(
            run_simulation,
            [(i, args, config) for i in range(1, args.num_simulations + 1)]
        )

    print("\nAll simulations completed successfully.")

# --------------------------------------------------
# Main CLI
# --------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Two-phase CNV-aware single-cell MDA simulation"
    )

    parser.add_argument("--num_simulations", type=int, default=1)
    parser.add_argument("--num_cores", type=int, default=multiprocessing.cpu_count())

    parser.add_argument("--config_file", required=True)
    parser.add_argument("--output_base", default="output")

    parser.add_argument("--patSeq_file", required=True)
    parser.add_argument("--matSeq_file", required=True)
    parser.add_argument("--cnv_bed_file", required=True)
    
    args = parser.parse_args()

    config = read_config(args.config_file, section="Simulation")
    print(f"Loaded configuration: {args.config_file}")

    run_all_simulations(args, config)

if __name__ == "__main__":
    main()
