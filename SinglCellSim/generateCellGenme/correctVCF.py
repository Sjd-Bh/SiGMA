import argparse
import os
import sys

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

def convert_to_vcf(input_vcf, output_vcf):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        # Write the VCF header
        outfile.write("##fileformat=VCFv4.2\n")
        outfile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
        outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        # Read input VCF-like data and convert to proper VCF format
        for line in infile:
            fields = line.strip().split()
            if len(fields) < 4:
                continue  # Skip lines that don't have enough columns

            chrom = fields[0]        # Chromosome or scaffold
            pos = fields[1]          # Position
            ref = fields[2]          # Reference allele
            alt = fields[3]          # Alternate allele
            
            # Format the line in proper VCF format
            vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tDP=.\tGT\t1/1\n"  # Placeholder for QUAL, FILTER, INFO, and FORMAT
            outfile.write(vcf_line)

def main():
    parser = argparse.ArgumentParser(description='Convert non-VCF formatted files into proper VCF format.')
    parser.add_argument('--files', nargs='+', help='List of VCF files to process', required=True)
    args = parser.parse_args()

    # Process each file
    for file_path in args.files:
        if os.path.isfile(file_path):
            # Get the base filename without the directory and extension
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            # Create a new output filename with "_edit.vcf"
            output_vcf = os.path.join(os.path.dirname(file_path), f"{base_name}_edit.vcf")
            print(f"Processing {file_path} -> {output_vcf}")
            convert_to_vcf(file_path, output_vcf)

if __name__ == "__main__":
    main()
