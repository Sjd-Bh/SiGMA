import argparse
import pickle
import random
import os
import sys

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)


def read_fasta(filename):
    with open(filename, 'r') as f:
        name, seq = None, []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))

def write_fasta(filename, name, sequence, line_length=60):
    with open(filename, 'w') as f:
        f.write(f">{name}\n")
        # Write the sequence in chunks of 'line_length'
        for i in range(0, len(sequence), line_length):
            f.write(sequence[i:i + line_length] + "\n")

def read_vcf(vcf_file):
    snvs = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            chrom, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
            snvs.append((chrom, pos, ref, alt))
    return snvs

def makeCell(pat_file, mat_file, pkl_file, snvs_file, select_cell, output_folder):
    # Read coalescent tree from pkl file
    with open(pkl_file, 'rb') as f:
        coalescent_tree = pickle.load(f)

    # Read SNVs from VCF file
    snvs = read_vcf(snvs_file)
    snvs = [(chrom, pos - 1, ref, alt) for (chrom, pos, ref, alt) in snvs]

    # Read paternal and maternal genomes
    pat_name, pat_seq = next(read_fasta(pat_file))
    mat_name, mat_seq = next(read_fasta(mat_file))

    # Introduce SNVs based on binomial distribution
    pat_snvs, mat_snvs = [], []
    selected_cell = f'cell_{select_cell}'
    cell_mutations = coalescent_tree['node_mutations'].get(selected_cell, [])
    for pos in cell_mutations:
        for snv in snvs:
            if snv[1] == pos:
                chrom, pos, ref, alt = snv
                if random.choice([True, False]):  # Randomly choose between paternal and maternal
                    pat_seq = pat_seq[:pos] + alt + pat_seq[pos:]
                    pat_snvs.append((chrom, pos + 1, ref, alt))
                else:
                    mat_seq = mat_seq[:pos] + alt + mat_seq[pos:]
                    mat_snvs.append((chrom, pos + 1, ref, alt))

    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Save updated genomes and SNVs to files
    write_fasta(os.path.join(output_folder, "paternal_cell.fasta"), pat_name, pat_seq)
    write_fasta(os.path.join(output_folder, "maternal_cell.fasta"), mat_name, mat_seq)

    with open(os.path.join(output_folder, "paternal_snvs.vcf"), 'w') as f:
        for snv in pat_snvs:
            f.write("\t".join(map(str, snv)) + "\n")
    with open(os.path.join(output_folder, "maternal_snvs.vcf"), 'w') as f:
        for snv in mat_snvs:
            f.write("\t".join(map(str, snv)) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--pat', required=True, help='Paternal genome FASTA file')
    parser.add_argument('--mat', required=True, help='Maternal genome FASTA file')
    parser.add_argument('--coal', required=True, help='Coalescent tree PKL file')
    parser.add_argument('--SNVs', required=True, help='VCF file with SNVs')
    parser.add_argument('--select-cell', type=int, required=True, help='Selected cell from the coalescent tree')
    parser.add_argument('--output', required=True, help='Output folder for generated files')
    args = parser.parse_args()

    makeCell(args.pat, args.mat, args.coal, args.SNVs, args.select_cell, args.output)

