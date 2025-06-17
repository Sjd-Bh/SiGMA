import pysam

def get_bulk_variants(bulk_vcf_path):
    bulk_vcf = pysam.VariantFile(bulk_vcf_path)
    bulk_positions = set()
    for rec in bulk_vcf:
        bulk_positions.add((rec.contig, rec.pos, str(rec.ref), str(rec.alts[0])))
    return bulk_positions

def check_single_cell_vafs(single_vcf_path, bulk_variants):
    vafs = {}
    keys = set()
    vcf_in = pysam.VariantFile(single_vcf_path)
    for rec in vcf_in:
        key = (rec.contig, rec.pos, str(rec.ref), str(rec.alts[0]))
        keys.add((rec.contig, rec.pos, str(rec.ref), str(rec.alts[0])))
        if key in bulk_variants:
            print("yes")
            sample = list(rec.samples.values())[0]
            if "AD" in sample and sample["AD"] and len(sample["AD"]) >= 2:
                ref_count, alt_count = sample["AD"][:2]
                total = ref_count + alt_count
                if total > 0:
                    vafs[key] = alt_count / total
                else:
                    vafs[key] = 0.0
    return vafs

import pickle

# Load the .pkl file
with open("../3nd_test_coal3/coalFolder/0.8Mb.pkl", "rb") as f:
    data = pickle.load(f)