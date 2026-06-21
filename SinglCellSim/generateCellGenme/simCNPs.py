import random
import argparse

def generate_het_cnps(fai_file, output_bed, num_cnps, min_len, max_len):
    # Read chromosome sizes from the fasta index file (.fai)
    chroms = {}
    with open(fai_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chroms[parts[0]] = int(parts[1])
    
    # Filter for standard chromosomes (e.g., chr1 to chr22, chrX) to avoid small contigs
    valid_chroms = [c for c in chroms.keys() if c.startswith('chr') and '_' not in c]
    if not valid_chroms:
        valid_chroms = list(chroms.keys()) # Fallback if no 'chr' prefix
        
    with open(output_bed, 'w') as out:
        for _ in range(num_cnps):
            chrom = random.choice(valid_chroms)
            
            # Ensure the CNP fits within the chromosome
            max_start = max(1, chroms[chrom] - max_len)
            start = random.randint(1, max_start)
            length = random.randint(min_len, max_len)
            end = start + length
            
            cnv_type = random.choice(['DEL', 'AMP'])
            allele = random.choice(['pat', 'mat']) # Assign to one allele to make it heterozygous
            
            # Write standard 5-column BED: chrom, start, end, type, allele
            out.write(f"{chrom}\t{start}\t{end}\t{cnv_type}\t{allele}\n")
            
    print(f"Generated {num_cnps} heterozygous CNPs in {output_bed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate random germline hetCNPs")
    parser.add_argument("--fai", required=True, help="Path to reference FASTA index (.fai)")
    parser.add_argument("--out", required=True, help="Output BED file path")
    parser.add_argument("--num", type=int, default=10, help="Number of CNPs to generate")
    parser.add_argument("--min_len", type=int, default=1000, help="Minimum CNP length")
    parser.add_argument("--max_len", type=int, default=50000, help="Maximum CNP length")
    args = parser.parse_args()
    
    generate_het_cnps(args.fai, args.out, args.num, args.min_len, args.max_len)
