#!/usr/bin/env python3
"""
baseline_vaf_caller.py
----------------------
A toy SNP caller: report biallelic SNVs whose VAF is in [min_vaf, max_vaf].

Requirements
------------
pip install pysam
Index the reference with `samtools faidx ref.fasta`
Index the BAM with `samtools index sample.bam`
"""

import argparse
import gzip
from collections import Counter
import pysam

VCF_HEADER = """##fileformat=VCFv4.2
##source=baseline_vaf_caller
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Ref and alt allele depths">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele fraction">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

def choose_alt_allele(base_counts, ref_base):
    """Return (alt_base, alt_depth) for the non-ref base with most support."""
    alt_base, depth = None, 0
    for base, cnt in base_counts.items():
        if base != ref_base.upper() and cnt > depth:
            alt_base, depth = base, cnt
    return alt_base, depth

def call_variants(bam_path, ref_path, min_vaf, max_vaf, min_depth, out_vcf):
    samfile  = pysam.AlignmentFile(bam_path, "rb")
    ref_fa   = pysam.FastaFile(ref_path)
    writer   = gzip.open(out_vcf, "wt") if out_vcf.endswith(".gz") else open(out_vcf, "w")
    writer.write(VCF_HEADER)

    for pileup_col in samfile.pileup(stepper="samtools", truncate=True):
        chrom = samfile.get_reference_name(pileup_col.reference_id)
        pos1  = pileup_col.reference_pos + 1            # 1-based
        ref   = ref_fa.fetch(chrom, pos1-1, pos1).upper()

        # Count A/C/G/T ignoring indels and low-quality bases
        bases = [p.alignment.query_sequence[p.query_position]
                 for p in pileup_col.pileups
                 if not p.is_del and not p.is_refskip]
        depth = len(bases)
        if depth < min_depth:
            continue

        base_counts = Counter(bases)
        alt, alt_depth = choose_alt_allele(base_counts, ref)
        if not alt:
            continue

        vaf = alt_depth / depth
        if min_vaf <= vaf <= max_vaf:
            info = f"DP={depth};AD={depth-alt_depth},{alt_depth};VAF={vaf:.3f}"
            writer.write(f"{chrom}\t{pos1}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")

    writer.close()
    samfile.close()
    ref_fa.close()
    print(f"Finished. Results written to {out_vcf}")

def parse_args():
    p = argparse.ArgumentParser(description="Baseline VAF-range SNP caller")
    p.add_argument("-b", "--bam", required=True, help="Input BAM/CRAM (indexed)")
    p.add_argument("-r", "--ref", required=True, help="Reference FASTA (indexed)")
    p.add_argument("--min-vaf", type=float, default=0.30, help="Lower VAF bound")
    p.add_argument("--max-vaf", type=float, default=0.70, help="Upper VAF bound")
    p.add_argument("--min-depth", type=int, default=8,  help="Min total depth")
    p.add_argument("-o", "--output", default="baseline_calls.vcf",
                   help="Output VCF (gz if .gz)")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    call_variants(args.bam, args.ref,
                  args.min_vaf, args.max_vaf,
                  args.min_depth, args.output)
