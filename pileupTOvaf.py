def calculate_vaf(pileup_file, output_file):
    with open(pileup_file, 'r') as file, open(output_file, 'w') as outfile:
        # Write the header to the output file
        outfile.write("Chromosome\tPosition\tReference\tDepth\tRefCount\tAltCount\tVAF\n")        
        for line in file:
            chrom, pos, ref, depth, bases, _ = line.strip().split('\t')[:6]
            ref_count = bases.count('.') + bases.count(',')
            alt_count = len(bases) - ref_count
            total_count = ref_count + alt_count
            vaf = alt_count / total_count if total_count > 0 else 0
            # Write the calculated values to the output file
            outfile.write(f"{chrom}\t{pos}\t{ref}\t{depth}\t{ref_count}\t{alt_count}\t{vaf:.2f}\n")
        print(f"VAF results saved to {output_file}")

# # Replace "chr10_snps.pileup" with your pileup file name and provide an output file name.
# calculate_vaf("chr10_snps.pileup", "chr10_snps_vaf_results.txt")