import argparse
import pickle
import random
import os
import subprocess
from collections import defaultdict
import numpy as np
from Bio import SeqIO

_BASES = ["A", "C", "G", "T"]

# ------------------------------
# Utility functions
# ------------------------------

def load_tree_pickle(pkl_file):
    with open(pkl_file, "rb") as f:
        return pickle.load(f)

def choose_alt_base(ref_base):
    if ref_base is None or ref_base.upper() not in _BASES:
        return random.choice(_BASES)
    return random.choice([b for b in _BASES if b != ref_base.upper()])

def load_reference_base_at_positions(fasta_path, chrom, positions):
    """Return dict of 0-based pos -> base"""
    if fasta_path is None:
        return {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if rec.id == chrom or rec.id.split()[0] == chrom:
            seq = str(rec.seq).upper()
            return {p: seq[p] if 0 <= p < len(seq) else "N" for p in positions}
    return {}

def compress_and_index_vcf(vcf_path):
    """Run bgzip and tabix on the generated VCF file."""
    subprocess.run(["bgzip", "-f", vcf_path], check=True)
    gz_path = f"{vcf_path}.gz"
    subprocess.run([#"conda", "run", "-n", "picard",
	"tabix", "-p", "vcf", gz_path], check=True)
    return gz_path

# ------------------------------
# SNV assignment
# ------------------------------

def assign_snvs_and_propagate_allele_aware(tree, branch_lengths, genome_length, mutation_rate, pat_prob=0.5):
    mutations = {}
    node_mutations = {node: defaultdict(set) for node in tree.nodes}
    mutation_positions = set()
    snv_alt_map = {}

    for (parent, child), length in branch_lengths.items():
        num_mutations = np.random.poisson(mutation_rate * length * genome_length)
        branch_map = {}
        for _ in range(num_mutations):
            attempt = 0
            while True:
                attempt += 1
                pos = np.random.randint(0, genome_length)
                if pos not in mutation_positions:
                    mutation_positions.add(pos)
                    r = random.random()
                    if r < pat_prob / 2: # Proportional to user pat_prob
                        allele_set = {"pat"}
                    elif r < pat_prob:
                        allele_set = {"mat"}
                    else:
                        allele_set = {"pat", "mat"}
                    branch_map.setdefault(pos, set()).update(allele_set)
                    break
                if attempt > 10000:
                    raise RuntimeError("Could not place SNV; increase genome_length or lower mutation_rate")
        mutations[(parent, child)] = branch_map

    def propagate(node):
        for child in tree.successors(node):
            for pos, alleles in node_mutations[node].items():
                node_mutations[child][pos].update(alleles)
            for pos, alleles in mutations.get((node, child), {}).items():
                node_mutations[child][pos].update(alleles)
            propagate(child)

    root_candidates = [n for n in tree.nodes if tree.in_degree(n) == 0]
    root = root_candidates[0] if len(root_candidates) == 1 else (max(tree.nodes) if tree.nodes else None)
    node_mutations[root] = defaultdict(set)
    propagate(root)

    leaf_nodes = [n for n in tree.nodes if tree.out_degree(n) == 0]
    vaf_any, vaf_pat, vaf_mat, allele_tags = {}, {}, {}, {}

    for pos in mutation_positions:
        mut_pat = sum(1 for leaf in leaf_nodes if "pat" in node_mutations[leaf].get(pos,set()))
        mut_mat = sum(1 for leaf in leaf_nodes if "mat" in node_mutations[leaf].get(pos,set()))
        total_leaves = len(leaf_nodes)
        vaf_pat[pos] = mut_pat / total_leaves if total_leaves else 0
        vaf_mat[pos] = mut_mat / total_leaves if total_leaves else 0
        vaf_any[pos] = (mut_pat + mut_mat) / (2 * total_leaves) if total_leaves else 0

        if mut_pat > 0 and mut_mat > 0: allele_tags[pos] = "both"
        elif mut_pat > 0: allele_tags[pos] = "pat"
        elif mut_mat > 0: allele_tags[pos] = "mat"
        else: allele_tags[pos] = "."

    return mutations, node_mutations, vaf_any, vaf_pat, vaf_mat, allele_tags, snv_alt_map, mutation_positions

# ------------------------------
# CNV assignment
# ------------------------------

def assign_cnvs_and_propagate(tree, branch_lengths, genome_length, cnv_rate, mean_cnv_length=10000, max_attempts=1000):
    cnv_mutations = {}
    node_cnvs = {node: [] for node in tree.nodes}
    cnv_counter = 0

    for (parent, child), length in branch_lengths.items():
        num_cnvs = np.random.poisson(cnv_rate * length * genome_length)
        branch_cnvs = []
        for _ in range(num_cnvs):
            attempts = 0
            while True:
                attempts += 1
                start = int(np.random.randint(0, genome_length))
                raw_len = int(max(1, np.random.exponential(scale=mean_cnv_length)))
                end = min(genome_length, start + raw_len)
                if end <= start: continue
                
                cnv_type = random.choice(["DEL", "DUP"])
                
                # Assign allele origin to CNV
                r = random.random()
                if r < 0.33: allele = "pat"
                elif r < 0.66: allele = "mat"
                else: allele = "both"

                cid = f"cnv_{cnv_counter}"
                cnv_counter += 1
                branch_cnvs.append({"id": cid, "start": start, "end": end, "type": cnv_type, "allele": allele})
                break
                if attempts > max_attempts: break
        cnv_mutations[(parent, child)] = branch_cnvs

    def propagate_cnvs(node):
        for child in tree.successors(node):
            node_cnvs[child].extend(node_cnvs[node])
            node_cnvs[child].extend(cnv_mutations.get((node, child), []))
            propagate_cnvs(child)

    root_candidates = [n for n in tree.nodes if tree.in_degree(n) == 0]
    root = root_candidates[0] if len(root_candidates) == 1 else (max(tree.nodes) if tree.nodes else None)
    node_cnvs[root] = []
    propagate_cnvs(root)

    leaf_nodes = [n for n in tree.nodes if tree.out_degree(n) == 0]
    id_to_cnv = {}
    for branch_list in cnv_mutations.values():
        for c in branch_list:
            id_to_cnv[c["id"]] = (c["start"], c["end"], c["type"], c["allele"])

    cnv_info = {}
    for cid, (s,e,t, a) in id_to_cnv.items():
        prevalence = sum(1 for leaf in leaf_nodes if any(c["id"] == cid for c in node_cnvs[leaf])) / len(leaf_nodes) if leaf_nodes else 0
        cnv_info[cid] = {"prevalence": prevalence, "start": s, "end": e, "type": t, "allele": a}

    return cnv_mutations, node_cnvs, cnv_info

# ------------------------------
# Per-Cell Output Generators
# ------------------------------

def write_per_cell_outputs(out_dir, leaf_nodes, node_mutations, node_cnvs, positions, snv_alt_map, paeref_bases, matref_bases, chrom="chr1"):
    os.makedirs(out_dir, exist_ok=True)
    
    for cell in leaf_nodes:
        # 1. Write Cell VCF
        cell_vcf_path = os.path.join(out_dir, f"{cell}_snvs.vcf")
        with open(cell_vcf_path, "w") as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write(f"##contig=<ID={chrom}>\n")
            vcf.write('##INFO=<ID=ALLELES,Number=.,Type=String,Description="Allele origins (pat,mat,both)">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{cell}\n")

            cell_muts = node_mutations.get(cell, {})
            for pos in sorted(cell_muts.keys()):
                alleles = cell_muts[pos]
                if not alleles: continue
                
                pos1 = pos + 1
                pat_ref = paeref_bases.get(pos, "N") if paeref_bases else "N"
                mat_ref = matref_bases.get(pos, "N") if matref_bases else "N"
                ref = pat_ref if paeref_bases else (mat_ref if matref_bases else "N")
                alt = snv_alt_map.get(pos)

                if "pat" in alleles and "mat" in alleles:
                    gt = "1|1"
                    allele_str = "both"
                else:
                    # Randomly assign phased het genotypes
                    gt = random.choice(["0|1", "1|0"])
                    allele_str = "pat" if "pat" in alleles else "mat"

                vcf.write(f"{chrom}\t{pos1}\t.\t{ref}\t{alt}\t.\tPASS\tALLELES={allele_str}\tGT\t{gt}\n")
        
        # Compress and index the cell VCF
        compress_and_index_vcf(cell_vcf_path)

        # 2. Write Cell CNV BED
        cell_bed_path = os.path.join(out_dir, f"{cell}_cnvs.bed")
        with open(cell_bed_path, "w") as bed:
            #bed.write("#chrom\tstart\tend\tid\ttype\tallele\n")
            for cnv in node_cnvs.get(cell, []):
                bed.write(f"{chrom}\t{cnv['start']}\t{cnv['end']}\t{cnv['id']}\t{cnv['type']}\t{cnv['allele']}\n")

# ------------------------------
# Main function
# ------------------------------

def main(args):
    data = load_tree_pickle(args.input)
    tree = data.get("tree")
    branch_lengths = data.get("branch_lengths") or {}
    if not branch_lengths:
        for u,v,d in tree.edges(data=True):
            branch_lengths[(u,v)] = float(d.get("time", 1.0))

    genome_length = args.genome_length or data.get("genome_length")

    # SNVs & CNVs
    mutations, node_mutations, vaf_any, vaf_pat, vaf_mat, allele_tags, snv_alt_map, mutation_positions = assign_snvs_and_propagate_allele_aware(
        tree, branch_lengths, genome_length, args.mutation_rate, pat_prob=args.pat_prob)

    cnv_mutations, node_cnvs, cnv_info = assign_cnvs_and_propagate(
        tree, branch_lengths, genome_length, args.cnv_rate, mean_cnv_length=args.mean_cnv_length)

    positions = sorted(list(mutation_positions))
    paeref_bases = load_reference_base_at_positions(args.patref, args.chrom, positions) if args.patref else {}
    matref_bases = load_reference_base_at_positions(args.matref, args.chrom, positions) if args.matref else {}

    for pos in positions:
        if pos not in snv_alt_map:
            pref_ref = paeref_bases.get(pos) if paeref_bases else (matref_bases.get(pos) if matref_bases else None)
            snv_alt_map[pos] = choose_alt_base(pref_ref)

    # Output to pickle
    # Output to pickle
    data_out = dict(data)
    data_out.update({
        "mutations": mutations,
        # Update this line to pack the alleles AND the alt_base into a tuple!
        "node_mutations": {n: {p: (set(a), snv_alt_map[p]) for p, a in d.items()} for n, d in node_mutations.items()},
        "cnv_mutations": cnv_mutations,
        "node_cnvs": node_cnvs,
        "cnv_info": cnv_info,
    })

    with open(args.output, "wb") as f:
        pickle.dump(data_out, f)
    print(f"[INFO] Updated pickle written to: {args.output}")

    # Generate Per-Cell VCFs and BEDs
    leaf_nodes = [n for n in tree.nodes if tree.out_degree(n) == 0]
    if args.cell_out_dir:
        write_per_cell_outputs(
            args.cell_out_dir, leaf_nodes, node_mutations, node_cnvs, 
            positions, snv_alt_map, paeref_bases, matref_bases, chrom=args.chrom
        )
        print(f"[INFO] Cell-specific VCFs (bgzipped/tabixed) and BEDs written to: {args.cell_out_dir}")

    print(f"[SUMMARY] SNVs: {len(positions)} | CNVs: {len(cnv_info)} | Leaves: {len(leaf_nodes)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--genome_length", type=int, required=True)
    parser.add_argument("--mutation_rate", type=float, default=1e-8)
    parser.add_argument("--cnv_rate", type=float, default=1e-9)
    parser.add_argument("--mean_cnv_length", type=int, default=10000)
    parser.add_argument("--cell_out_dir", type=str, default="cell_outputs", help="Directory for per-cell VCFs and BEDs")
    parser.add_argument("--patref", type=str, default=None)
    parser.add_argument("--matref", type=str, default=None)
    parser.add_argument("--chrom", type=str, default="chr1")
    parser.add_argument("--pat_prob", type=float, default=0.5)
    args = parser.parse_args()
    main(args)
