###########################################
import argparse
import os
from configs.configFunctions import read_config, read_fasta
from singleCellGenomeSim.singleCellGenomeGen import generate_single_cell_genome_from_fasta, save_to_vcf,save_cnvs_to_bed, save_fasta


if os.name == 'nt':  # Windows
    os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    params = read_config("configs\\generateCellParam.ini", "Parameters")
    templates = read_fasta('AmpSim\\templates.fasta')
    output_dir = "output\\"  # Define your output directory for Windows
    os.makedirs(output_dir, exist_ok=True)  # Create the output directory if it doesn't exist
    # Generate single-cell genomes
    paternal_genome, snps_paternal, snvs_paternal, cnvs_paternal = generate_single_cell_genome_from_fasta(
        templates, genome_name="paternal", num_snps=int(params["num_snp"]), num_snvs=int(params["num_snv"]), num_cnvs=int(params["num_cnv"]))
    maternal_genome, snps_maternal, snvs_maternal, cnvs_maternal = generate_single_cell_genome_from_fasta(
        templates, genome_name="maternal", num_snps=int(params["num_snp"]), num_snvs=int(params["num_snv"]), num_cnvs=int(params["num_cnv"]))

    save_fasta(paternal_genome,"Paternal_Genome", "paternal_genome.fasta")
    save_fasta(maternal_genome,"Maternal_Genome", "maternal_genome.fasta")
else:  # Linux
    parser = argparse.ArgumentParser(description='Generate Single Cell Genome')
    parser.add_argument('--ref', help='Path to the reference genome file (reference.fasta)')
    parser.add_argument('--num-snp', type=int, help='Number of SNPs to introduce')
    parser.add_argument('--num-snv', type=int, help='Number of SNVs to introduce')
    parser.add_argument('--num-cnv', type=int, help='Number of CNVs to introduce')
    parser.add_argument('--o', dest='output', help='Output folder for results')
    
    args = parser.parse_args()

    # Create the output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)

    # Ensure output directory ends with a trailing slash
    if not args.output.endswith('/'):
        args.output += '/'

    #params = read_config("configs/config.ini", "Parameters")
    templates = read_fasta(args.ref)
    output_dir = args.output

    # Generate single-cell genomes
    paternal_genome, snps_paternal, snvs_paternal, cnvs_paternal = generate_single_cell_genome_from_fasta(
        templates,genome_name = "paternal" , num_snps=args.num_snp, num_snvs=args.num_snv, num_cnvs=args.num_cnv)
    maternal_genome, snps_maternal, snvs_maternal, cnvs_maternal = generate_single_cell_genome_from_fasta(
        templates, genome_name = "maternal" ,num_snps=args.num_snp, num_snvs=args.num_snv, num_cnvs=args.num_cnv)

    # Combine paternal and maternal genomes into a single-cell genome (diploid)
    # single_cell_genome = ">Single_Cell_Genome\n" + paternal_genome + maternal_genome

    # Save the combined single-cell genome as a FASTA file
    save_fasta(paternal_genome, "Paternal_Genome", os.path.join(output_dir, 'paternal_genome.fasta'))
    save_fasta( maternal_genome,"Maternal_Genome", os.path.join(output_dir, 'maternal_genome.fasta'))


    # Save SNPs, SNVs, and CNVs to separate files
    save_to_vcf(snps_paternal, os.path.join(output_dir, 'snps_paternal.vcf'), paternal_genome)
    save_to_vcf(snps_maternal, os.path.join(output_dir, 'snps_maternal.vcf'), maternal_genome)
    save_to_vcf(snvs_paternal, os.path.join(output_dir, 'snvs_paternal.vcf'), paternal_genome)
    save_to_vcf(snvs_maternal, os.path.join(output_dir, 'snvs_maternal.vcf'), maternal_genome)
    save_cnvs_to_bed(cnvs_paternal, os.path.join(output_dir, 'cnvs_paternal.bed'), paternal_genome)
    save_cnvs_to_bed(cnvs_maternal, os.path.join(output_dir, 'cnvs_maternal.bed'), maternal_genome)

    print("Process completed successfully!")

  












