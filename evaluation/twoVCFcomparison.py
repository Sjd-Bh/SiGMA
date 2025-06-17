# twoVCFcomparison.py

import argparse
import csv
import gzip
import os

def open_vcf(vcf_file):
    if vcf_file.endswith(".gz"):
        return gzip.open(vcf_file, 'rt')
    else:
        return open(vcf_file, 'r')

def read_vcf_positions(vcf_file):
    positions = set()
    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])  # VCF is 1-based
            ref = fields[3]
            alt = fields[4]
            positions.add((chrom, pos, ref, alt))
    return positions

def write_tp_vcf(called_vcf_file, tp_set, output_vcf):
    with open_vcf(called_vcf_file) as f_in, open(output_vcf, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
            else:
                fields = line.strip().split('\t')
                key = (fields[0], int(fields[1]), fields[3], fields[4])
                if key in tp_set:
                    f_out.write(line)

def compare_variants(ground_truth, called_variants, tp_output_vcf):
    gt_set = read_vcf_positions(ground_truth)
    called_set = read_vcf_positions(called_variants)

    tp_set = gt_set & called_set
    fp = len(called_set - gt_set)
    fn = len(gt_set - called_set)

    write_tp_vcf(called_variants, tp_set, tp_output_vcf)
    return len(tp_set), fp, fn

def main():
    parser = argparse.ArgumentParser(description='Compare called variants to ground truth and compute TP, FP, FN.')
    parser.add_argument('--ground-truth', required=True, help='Ground truth VCF file')
    parser.add_argument('--called-var', required=True, help='Called variants VCF file')
    parser.add_argument('--output', required=True, help='Output CSV file for TP/FP/FN counts')
    parser.add_argument('--tp-vcf', required=True, help='Output VCF file to save TP variants')

    args = parser.parse_args()

    tp, fp, fn = compare_variants(args.ground_truth, args.called_var, args.tp_vcf)

    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['TP', 'FP', 'FN'])
        writer.writeheader()
        writer.writerow({'TP': tp, 'FP': fp, 'FN': fn})

    print(f"TP, FP, FN counts saved to {args.output}")
    print(f"TP variants saved to {args.tp_vcf}")

if __name__ == '__main__':
    main()
