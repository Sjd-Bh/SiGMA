import argparse
import os
import sys
import subprocess

# Ensure the script is run from the project root directory
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))
sys.path.insert(0, project_root)

def convert_to_vcf(input_vcf, output_vcf):
    with open(input_vcf, 'r') as infile:
        # Read the first line to get the contig name
        first_line = infile.readline().strip()
        fields = first_line.split()
        if len(fields) < 4:
            raise ValueError(f"Input file {input_vcf} does not have enough columns in the first line.")
        
        contig = fields[0]  # Extract the contig name (e.g., '400kb')

        # Write the VCF header
        with open(output_vcf, 'w') as outfile:
            outfile.write("##fileformat=VCFv4.2\n")
            outfile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
            outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            outfile.write(f"##contig=<ID={contig}>\n")  # Use the extracted contig name
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

            # Write the rest of the data
            outfile.write(first_line + '\n')  # Write the first line again
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

def process_vcf_file(file_path):
    # Get the base filename without the directory and extension
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    # Create a new output filename with "_edit.vcf"
    output_vcf = os.path.join(os.path.dirname(file_path), f"{base_name}_edit.vcf")
    print(f"Processing {file_path} -> {output_vcf}")

    convert_to_vcf(file_path, output_vcf)

    # Sort the VCF file
    sorted_vcf = os.path.join(os.path.dirname(file_path), f"{base_name}_sorted.vcf")
    print(f"Sorting {output_vcf} -> {sorted_vcf}")
    subprocess.run(['bcftools', 'sort', '-o', sorted_vcf, output_vcf], check=True)

    # Compress the sorted VCF file
    compressed_vcf = f"{sorted_vcf}.gz"
    print(f"Compressing {sorted_vcf} -> {compressed_vcf}")
    subprocess.run(['bgzip', sorted_vcf], check=True)

    # Index the compressed VCF file
    print(f"Indexing {compressed_vcf}")
    subprocess.run(['bcftools', 'index', compressed_vcf], check=True)

def main():
    parser = argparse.ArgumentParser(description='Convert non-VCF formatted files into proper VCF format.')
    parser.add_argument('--files', nargs='+', help='List of VCF files to process', required=True)
    args = parser.parse_args()

    # Process each file
    for file_path in args.files:
        if os.path.isfile(file_path):
            process_vcf_file(file_path)
        else:
            print(f"File not found: {file_path}")

if __name__ == "__main__":
    main()



# input_vcf = '../../test\\out.vcf'  # Your non-standard VCF file
# output_vcf = '../../test\\out_edit.vcf'  # The correct VCF file to be created
# convert_to_vcf(input_vcf, output_vcf)