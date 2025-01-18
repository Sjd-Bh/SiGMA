def merge_vcf(file1, file2, output_file):
    def parse_vcf(filename):
        """Parse a VCF file and separate header and variants."""
        header = []
        variants = []
        with open(filename, "r") as f:
            for line in f:
                if line.startswith("#"):
                    header.append(line.strip())
                else:
                    variants.append(line.strip())
        return header, variants
    # Parse both VCF files
    header1, variants1 = parse_vcf(file1)
    header2, variants2 = parse_vcf(file2)
    # Ensure headers match
    if header1 != header2:
        raise ValueError("Headers of the two VCF files do not match.")
    # Combine and sort variants
    all_variants = list(set(variants1 + variants2))  # Remove duplicates
    all_variants.sort(key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))  # Sort by CHROM, POS
    # Write merged VCF
    with open(output_file, "w") as out:
        out.write("\n".join(header1) + "\n")  # Write header
        out.write("\n".join(all_variants) + "\n")  # Write sorted variants
