# calculate_ADOs.py

import argparse
import os
import sys
import subprocess
import tempfile
import csv

import twoVCFcomparison 

def extract_homozygous(tp_vcf, temp_dir):
    hom_vcf_gz = tp_vcf.replace('.vcf', '_hom.vcf.gz')
    subprocess.run(['bcftools', 'view', '-g', 'hom', '-Oz', '-o', hom_vcf_gz, tp_vcf], check=True)
    subprocess.run(['bcftools', 'index', '-f', hom_vcf_gz], check=True)

    count = int(subprocess.check_output(
        f'bcftools view {hom_vcf_gz} | grep -v "^#" | wc -l',
        shell=True
    ).decode().strip())

    return count

def main():
    parser = argparse.ArgumentParser(description='Calculate ADOs from multiple VCF files.')
    parser.add_argument('--truth-vcf', required=True, help='Ground truth VCF file')
    parser.add_argument('--called-vcfs', nargs='+', required=True, help='List of called VCF files (expanded by shell)')
    parser.add_argument('--output', required=True, help='Output CSV file for ADO summary')

    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as temp_dir, open(args.output, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['File', 'FN', 'Homozygous_TP', 'ADO'])
        writer.writeheader()

        for vcf_file in sorted(args.called_vcfs):
            try:
                base = os.path.splitext(os.path.basename(vcf_file))[0].replace('.vcf', '').replace('.gz', '')
                tp_vcf = os.path.join(temp_dir, f'{base}_TP.vcf')

                tp, fp, fn = twoVCFcomparison.compare_variants(args.truth_vcf, vcf_file, tp_vcf)
                hom_tp_count = extract_homozygous(tp_vcf, temp_dir)
                ado = fn + hom_tp_count

                writer.writerow({'File': os.path.basename(vcf_file), 'FN': fn, 'Homozygous_TP': hom_tp_count, 'ADO': ado})
                print(f"{vcf_file} → FN={fn}, Homo_TP={hom_tp_count}, ADO={ado}")
            except Exception as e:
                print(f"Error processing {vcf_file}: {e}")

if __name__ == '__main__':
    main()
