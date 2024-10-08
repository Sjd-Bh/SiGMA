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

        # Check if the first line has enough columns to extract the contig name
        if len(fields) < 4:
            raise ValueError(f"Input file {input_vcf} does not have enough columns in the first line.")
        
        # Extract the contig name from the first line
        contig = fields[0]  # Use the first field as the contig name

        # Write the VCF header
        with open(output_vcf, 'w') as outfile:
            outfile.write("##fileformat=VCFv4.2\n")
            outfile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
            outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            outfile.write(f"##contig=<ID={contig}>\n")  # Add the contig name to the VCF header
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

            # Write the first line as is, but ensure it has enough fields
            if len(fields) >= 4:
                chrom, pos, ref, alt = fields[0], fields[1], fields[2], fields[3]
                vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tDP=.\tGT\t1/1\n"
                outfile.write(vcf_line)

            # Process the rest of the lines
            for line in infile:
                fields = line.strip().split()
                if len(fields) < 4:
                    print(f"Skipping line with insufficient fields: {line.strip()}")
                    continue  # Skip lines that don't have enough columns

                # Extract fields safely
                chrom = fields[0] if len(fields) > 0 else "."
                pos = fields[1] if len(fields) > 1 else "."
                ref = fields[2] if len(fields) > 2 else "."
                alt = fields[3] if len(fields) > 3 else "."
                
                # Format the line in proper VCF format
                vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tDP=.\tGT\t1/1\n"
                outfile.write(vcf_line)

def process_vcf_file(file_path):
    # Get the base filename without the directory and extension
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    # Create a new output filename with "_edit.vcf"
    output_vcf = os.path.join(os.path.dirname(file_path), f"{base_name}_edit.vcf")
    print(f"Processing {file_path} -> {output_vcf}")

    # Convert the file
    convert_to_vcf(file_path, output_vcf)
    
    # Verify if the VCF file was created
    if not os.path.exists(output_vcf):
        print(f"Error: {output_vcf} not created.")
        return

    # Sort the VCF file
    sorted_vcf = os.path.join(os.path.dirname(file_path), f"{base_name}_sorted.vcf")
    print(f"Sorting {output_vcf} -> {sorted_vcf}")
    try:
        subprocess.run(['bcftools', 'sort', '-o', sorted_vcf, output_vcf], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error sorting VCF file: {e}")
        return

    # Compress the sorted VCF file
    compressed_vcf = f"{sorted_vcf}.gz"
    print(f"Compressing {sorted_vcf} -> {compressed_vcf}")
    try:
        subprocess.run(['bgzip', sorted_vcf], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error compressing VCF file: {e}")
        return

    # Index the compressed VCF file
    print(f"Indexing {compressed_vcf}")
    try:
        subprocess.run(['bcftools', 'index', compressed_vcf], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error indexing VCF file: {e}")
        return

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