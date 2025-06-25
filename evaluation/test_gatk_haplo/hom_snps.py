import subprocess
import os
from collections import defaultdict
import argparse

# def run_gatk_bulk_haplotypecaller(bam, reference, output_vcf):
#     cmd = [
#         "gatk", "HaplotypeCaller",
#         "-R", reference,
#         "-I", bam,
#         "-O", output_vcf,
#         "--standard-min-confidence-threshold-for-calling", "30.0",
#         "--dont-use-soft-clipped-bases",
#         "--native-pair-hmm-threads","5"
#     ]
#     subprocess.run(cmd, check=True)
#     print(f"Bulk variant calling completed: {output_vcf}")

def run_gatk_sc_haplotypecaller(bam, reference, output_vcf):
    cmd = [
        "gatk", "HaplotypeCaller",
        "-R", reference,
        "-I", bam,
        "-O", output_vcf,
        "--standard-min-confidence-threshold-for-calling", "1.0",
        "--pcr-indel-model", "NONE",
        "--dont-use-soft-clipped-bases",
        "--native-pair-hmm-threads","5"
    ]
    subprocess.run(cmd, check=True)
    print(f"Single-cell variant calling completed: {output_vcf}")

# def extract_homozygous_snps(input_vcf, output_vcf):
#     shell_cmd = f"""
#     conda activate bcftools && \
#     bcftools view -g hom -v snps -o {output_vcf} -O v {input_vcf}
#     """

#     # Run the command in shell mode
#     subprocess.run(shell_cmd, shell=True, executable='/bin/bash', check=True)
#     print(f"Homozygous SNPs extracted to: {output_vcf}")

def load_genotypes(vcf_file):
    genotypes = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            chrom, pos, ref, alt, fmt, sample_data = fields[0], fields[1], fields[3], fields[4], fields[8], fields[9]
            fmt_keys = fmt.split(":")
            sample_vals = sample_data.split(":")
            fmt_dict = dict(zip(fmt_keys, sample_vals))
            gt = fmt_dict.get("GT")
            if gt:
                alleles = gt.replace('|', '/').split('/')
                genotypes[(chrom, pos)] = alleles
    return genotypes

def calculate_misclassification_rate(bulk_hom_vcf, sc_vcf, output_txt):
    bulk_snps = load_genotypes(bulk_hom_vcf)
    sc_snps = load_genotypes(sc_vcf)

    total_hom_snps = len(bulk_snps)
    miscalled_as_het = 0

    for key, gt in bulk_snps.items():
        sc_gt = sc_snps.get(key)
        if sc_gt and sc_gt[0] != sc_gt[1]:  # heterozygous in single cell
            miscalled_as_het += 1

    percent = (miscalled_as_het / total_hom_snps) * 100 if total_hom_snps else 0
    result = (
        f"Misclassified homozygous SNPs as heterozygous: "
        f"{miscalled_as_het}/{total_hom_snps} ({percent:.2f}%)"
    )

    with open(output_txt, "w") as out:
        out.write(result + "\n")

    print(result)

def main():
    parser = argparse.ArgumentParser(description="GATK variant calling and misclassification evaluation")
    parser.add_argument('--bulk-hom-vcf', required=True)
    parser.add_argument('--sc-bam', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--prefix', required=True, help="Base name for output files (e.g., sample1)")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # bulk_vcf = os.path.join(args.outdir, f"{args.prefix}_bulk.vcf")
    # bulk_hom_vcf = os.path.join(args.outdir, f"{args.prefix}_bulk_hom.vcf")
    sc_vcf = os.path.join(args.outdir, f"{args.prefix}_sc.vcf")
    result_txt = os.path.join(args.outdir, f"{args.prefix}_misclassification.txt")

    # run_gatk_bulk_haplotypecaller(args.bulk_bam, args.reference, bulk_vcf)
    # extract_homozygous_snps(bulk_vcf, bulk_hom_vcf)
    run_gatk_sc_haplotypecaller(args.sc_bam, args.reference, sc_vcf)

    calculate_misclassification_rate(args.bulk_hom_vcf, sc_vcf, result_txt)

if __name__ == "__main__":
    main()
