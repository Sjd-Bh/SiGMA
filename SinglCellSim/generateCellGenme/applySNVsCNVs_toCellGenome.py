#!/usr/bin/env python3
import argparse
import pickle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

# -----------------------------
# Utility functions
# -----------------------------

def load_fasta_as_dict(fasta_file):
    """Load a FASTA file into a dictionary {chrom: sequence}."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = list(str(record.seq))  # list of chars for easy mutation
    return sequences


def save_fasta_from_dict(sequences, output_file, prefix=""):
    """Save a dictionary {chrom: [bases]} to a FASTA file."""
    records = []
    for chrom, seq_list in sequences.items():
        seq_str = "".join(seq_list)
        records.append(SeqRecord(Seq(seq_str), id=f"{chrom}", description=""))
    SeqIO.write(records, output_file, "fasta")
    print(f"[INFO] Saved FASTA: {output_file}")


def mutate_sequence(sequence_list, position, new_base):
    """Mutate a sequence (list) at a 1-based position."""
    pos_index = position - 1
    if 0 <= pos_index < len(sequence_list):
        sequence_list[pos_index] = new_base.upper()
    return sequence_list


def apply_snvs(sequence_dict, snvs, genome_type=None):
    """Apply SNVs to genome dictionary."""
    count = 0
    chrom_default = list(sequence_dict.keys())[0]

    for snv in snvs:
        if isinstance(snv, int):
            mutate_sequence(sequence_dict[chrom_default], snv + 1, "A")
            count += 1
            continue

        targets = snv.get("target", [])
        if genome_type and genome_type not in targets and "both" not in targets:
            continue

        chrom = snv.get("chrom", chrom_default)
        if chrom not in sequence_dict:
            continue

        pos = snv["pos"] + 1
        alt_base = snv.get("alt", "A")
        mutate_sequence(sequence_dict[chrom], pos, alt_base)
        count += 1

    print(f"[INFO] Applied {count} SNVs to {genome_type} genome")
    return sequence_dict


def apply_cnvs(sequence_dict, cnvs, genome_type=None):
    """Apply CNVs (duplications or deletions) to genome."""
    total_cnvs = 0
    chrom_default = list(sequence_dict.keys())[0]

    for cnv in cnvs:
        chrom = cnv.get("chrom", chrom_default)
        target = cnv.get("target", "both")
        cnv_type = cnv.get("type", "").lower()

        # Normalize CNV type abbreviations
        if cnv_type in ["dup", "duplication"]:
            cnv_type = "duplication"
        elif cnv_type in ["del", "deletion"]:
            cnv_type = "deletion"
        else:
            continue  # skip unknown type

        if genome_type and target not in (genome_type, "both"):
            continue
        if chrom not in sequence_dict:
            continue

        start = int(cnv["start"])
        end = min(int(cnv["end"]), len(sequence_dict[chrom]))

        if cnv_type == "deletion":
            for i in range(start, end):
                sequence_dict[chrom][i] = "N"
        elif cnv_type == "duplication":
            segment = sequence_dict[chrom][start:end]
            sequence_dict[chrom] = sequence_dict[chrom][:end] + segment + sequence_dict[chrom][end:]

        total_cnvs += 1

    print(f"[INFO] Applied {total_cnvs} CNVs to {genome_type} genome")
    return sequence_dict


def merge_genomes(pat_dict, mat_dict):
    """Merge paternal and maternal sequences into a single diploid genome."""
    merged = {}
    for chrom in pat_dict:
        if chrom in mat_dict:
            merged_seq = []
            for p, m in zip(pat_dict[chrom], mat_dict[chrom]):
                merged_seq.append(p if p == m else f"{p}|{m}")
            merged[chrom] = merged_seq
    return merged


# -----------------------------
# Main
# -----------------------------

def generate_single_cell(pat_fasta, mat_fasta, pkl_file, node, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    node = int(node)

    print("[INFO] Loading genomes...")
    pat_genome = load_fasta_as_dict(pat_fasta)
    mat_genome = load_fasta_as_dict(mat_fasta)

    print("[INFO] Loading coalescent PKL...")
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    if node not in data["node_mutations"]:
        raise ValueError(f"Node '{node}' not found in PKL file")

    node_snvs = data["node_mutations"][node]
    # You need to adjust this based on what your pickle file actually contains!
    snvs = [{"pos": pos, "target": list(targets), "alt": alt_base} for pos, (targets, alt_base) in node_snvs.items()]
    node_cnvs = data.get("node_cnvs", {}).get(node, [])
    print(f"[DEBUG] CNVs for node {node}: {node_cnvs}")

    print(f"[INFO] Applying SNVs and CNVs to node {node}...")

    apply_snvs(pat_genome, snvs, genome_type="pat")
    apply_snvs(mat_genome, snvs, genome_type="mat")

    #apply_cnvs(pat_genome, node_cnvs, genome_type="pat")
    #apply_cnvs(mat_genome, node_cnvs, genome_type="mat")

    pat_output = os.path.join(output_dir, f"{node}_paternal.fasta")
    mat_output = os.path.join(output_dir, f"{node}_maternal.fasta")
    merged_output = os.path.join(output_dir, f"{node}_merged.fasta")

    save_fasta_from_dict(pat_genome, pat_output, prefix="pat_")
    save_fasta_from_dict(mat_genome, mat_output, prefix="mat_")

    merged_genome = merge_genomes(pat_genome, mat_genome)
    save_fasta_from_dict(merged_genome, merged_output, prefix="merged_")

    print(f"[INFO] Single-cell genome generation completed for node {node}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate single-cell genome with SNVs/CNVs applied")
    parser.add_argument("--pat", required=True)
    parser.add_argument("--mat", required=True)
    parser.add_argument("--pkl", required=True)
    parser.add_argument("--node", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    generate_single_cell(args.pat, args.mat, args.pkl, args.node, args.output)
