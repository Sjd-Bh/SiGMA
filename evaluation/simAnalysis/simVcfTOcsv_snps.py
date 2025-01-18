import pysam
import csv

def calculate_vaf(bam_file, vcf_file, output_file):
    def parse_vcf(vcf_filename):
        """Parse the VCF file and return a list of SNPs."""
        snps = []
        with open(vcf_filename, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    snps.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt})
        return snps
    def count_alleles(bam, chrom, pos, ref, alt):
        """Count REF and ALT alleles at a given position."""
        ref_count = 0
        alt_count = 0
        total_count = 0        
        # Fetch reads overlapping the position
        for pileupcolumn in bam.pileup(chrom, pos - 1, pos, truncate=True):
            if pileupcolumn.pos == pos - 1:  # Match position (0-based in pysam)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:  # Exclude deletions/skips
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if base == ref:
                            ref_count += 1
                        elif base == alt:
                            alt_count += 1
                        total_count += 1
        return ref_count, alt_count, total_count
    # Parse VCF file
    snps = parse_vcf(vcf_file)
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    # Calculate VAF for each SNP
    results = []
    for snp in snps:
        ref_count, alt_count, total_count = count_alleles(bam, snp["chrom"], snp["pos"], snp["ref"], snp["alt"])
        vaf = alt_count / total_count if total_count > 0 else 0
        results.append({"chrom": snp["chrom"], "pos": snp["pos"], "ref": snp["ref"], "alt": snp["alt"], 
                        "ref_count": ref_count, "alt_count": alt_count, "vaf": vaf})
    bam.close()
    # Write results to CSV
    with open(output_file, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["CHROM", "POS", "REF", "ALT", "REF_COUNT", "ALT_COUNT", "VAF"])  # Header
        for result in results:
            csvwriter.writerow([
                result['chrom'], result['pos'], result['ref'], result['alt'], 
                result['ref_count'], result['alt_count'], f"{result['vaf']:.4f}"
            ])

# Example usage
# bam_file = "single_cell.bam"
# vcf_file = "merged.vcf"
# output_file = "vaf_results.csv"

# calculate_vaf(bam_file, vcf_file, output_file)
