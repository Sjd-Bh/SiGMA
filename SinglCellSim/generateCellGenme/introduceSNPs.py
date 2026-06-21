import argparse
from Bio import SeqIO

def read_fasta_as_bytearray(fasta_path):
    record = SeqIO.read(fasta_path, "fasta")
    return bytearray(str(record.seq).upper().encode())

def read_vcf_file(vcf_path):
    mutations = {}
    # Use gzip or standard open depending on file extension
    open_func = open
    if str(vcf_path).endswith('.gz'):
        import gzip
        open_func = gzip.open

    with open_func(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            
            ref = cols[3]
            alt = cols[4]
            
            # Safety checks: only process standard single-base substitutions
            if len(ref) != 1 or len(alt) != 1 or ',' in alt:
                continue
                
            pos = int(cols[1]) - 1  # convert to 0-based
            mutations[pos] = alt.encode()
    return mutations

def apply_snp_mutations_inplace(genome, mutations):
    for pos, alt in mutations.items():
        if pos < len(genome): # Ensure we don't index out of bounds
            genome[pos] = alt[0]   # SNP = single base
    return genome

def write_fasta(genome, output_path, header="mutated_genome", line_width=70):
    with open(output_path, "wb") as f:
        f.write(f">{header}\n".encode())
        for i in range(0, len(genome), line_width):
            f.write(genome[i:i + line_width] + b"\n")

def process_genome_with_snp_mutations(fasta_path, vcf_path, output_path, chrom):
    genome = read_fasta_as_bytearray(fasta_path)
    mutations = read_vcf_file(vcf_path)
    apply_snp_mutations_inplace(genome, mutations)
    write_fasta(genome, output_path, header=chrom)

def main():
    parser = argparse.ArgumentParser(
        description="Apply SNPs from VCF to a reference genome efficiently"
    )
    parser.add_argument("--ref", required=True)
    parser.add_argument("--vcf", required=True, help="Can be .vcf or .vcf.gz")
    parser.add_argument("--out", required=True)
    parser.add_argument("--chrom", required=True, help="Chromosome name for the FASTA header")

    args = parser.parse_args()

    process_genome_with_snp_mutations(args.ref, args.vcf, args.out, args.chrom)

if __name__ == "__main__":
    main()
