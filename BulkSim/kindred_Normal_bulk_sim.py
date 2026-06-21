import argparse
import os
import subprocess
import glob

############################################################
# Helper Functions
############################################################

def apply_bcftools_consensus(reference_fasta, vcf_file, output_fasta, haplotype=None):
    """Apply variants from a VCF file to a reference FASTA, fixing REF mismatches."""
    temp_vcf_gz = output_fasta + ".temp.vcf.gz"
    
    norm_cmd = [
        "conda", "run", "-n", "picard", "bcftools", "norm", 
        "--check-ref", "s", "-f", reference_fasta, 
        "-O", "z", "-o", temp_vcf_gz, vcf_file
    ]
    subprocess.run(norm_cmd, check=True, stderr=subprocess.DEVNULL)
    
    index_cmd = ["conda", "run", "-n", "picard", "bcftools", "index", "-t", temp_vcf_gz]
    subprocess.run(index_cmd, check=True, stderr=subprocess.DEVNULL)
    
    cons_cmd = ["conda", "run", "-n", "picard", "bcftools", "consensus", "-f", reference_fasta, temp_vcf_gz]
    if haplotype is not None:
        cons_cmd.extend(["-H", str(haplotype)])
    
    with open(output_fasta, "w") as out:
        subprocess.run(cons_cmd, stdout=out, check=True, stderr=subprocess.DEVNULL)
        
    if os.path.exists(temp_vcf_gz): os.remove(temp_vcf_gz)
    if os.path.exists(temp_vcf_gz + ".tbi"): os.remove(temp_vcf_gz + ".tbi")

def apply_cnvs_to_fasta(fasta_in, fasta_out, bed_files, target_allele):
    """Applies CNVs (DEL/AMP) from multiple BED files to FASTA based on target allele."""
    valid_beds = [b for b in bed_files if b and os.path.exists(b) and os.path.getsize(b) > 0]
    
    if not valid_beds:
        subprocess.run(["cp", fasta_in, fasta_out], check=True)
        return

    with open(fasta_in, 'r') as f:
        header = f.readline().strip()
        seq = "".join(line.strip() for line in f)
    
    cnvs = []
    for bed_file in valid_beds:
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue
                
                parts = line.split('\t')
                if len(parts) < 5: continue
                
                try:
                    start, end = int(parts[1]), int(parts[2])
                except ValueError:
                    continue 
                
                cnv_type, allele = parts[3], parts[4] 
                if allele == target_allele or allele == "both":
                    cnvs.append((start, end, cnv_type))
    
    if not cnvs:
        subprocess.run(["cp", fasta_in, fasta_out], check=True)
        return
                
    cnvs.sort(key=lambda x: x[0], reverse=True)
    for start, end, cnv_type in cnvs:
        if cnv_type == "DEL":
            seq = seq[:start] + seq[end:]
        elif cnv_type == "AMP":
            seq = seq[:start] + seq[start:end] + seq[start:]
            
    with open(fasta_out, 'w') as f:
        f.write(header + "\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

def run_art(fasta_file, coverage, read_length, insert_size, std_dev, output_prefix):
    """Run ART Illumina to simulate reads."""
    cmd = [
        "art_illumina", "-ss", "HS25", "-i", fasta_file, "-p",
        "-l", str(read_length), "-f", str(coverage),
        "-m", str(insert_size), "-s", str(std_dev), "-o", output_prefix
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)

def merge_fastqs(output_dir, bulk_prefix, clone_prefixes):
    """Merge multiple FASTQ files into one bulk FASTQ for R1 and R2."""
    with open(os.path.join(output_dir, f"{bulk_prefix}_R1.fq"), "wb") as out1:
        for prefix in clone_prefixes:
            if os.path.exists(f"{prefix}1.fq"):
                with open(f"{prefix}1.fq", "rb") as f: out1.write(f.read())

    with open(os.path.join(output_dir, f"{bulk_prefix}_R2.fq"), "wb") as out2:
        for prefix in clone_prefixes:
            if os.path.exists(f"{prefix}2.fq"):
                with open(f"{prefix}2.fq", "rb") as f: out2.write(f.read())

############################################################
# Main Simulation Workflow
############################################################

def main():
    parser = argparse.ArgumentParser(description="Simulate bulk and clone sequencing from per-cell VCFs and BEDs.")
    parser.add_argument("--paternal_fasta", required=True, help="Pre-generated baseline paternal FASTA")
    parser.add_argument("--maternal_fasta", required=True, help="Pre-generated baseline maternal FASTA")
    parser.add_argument("--het_cnp_bed", required=False, help="Optional: BED file containing germline/heterozygous CNPs")
    parser.add_argument("--cell_dir", required=True, help="Directory containing per-cell *_snvs.vcf.gz and *_cnvs.bed")
    parser.add_argument("--output_dir", required=True, help="Directory to save outputs")
    parser.add_argument("--tumor_coverage", type=float, default=30.0, help="Total tumor bulk AND individual clone coverage")
    parser.add_argument("--normal_coverage", type=float, default=30.0, help="Total normal bulk coverage")
    parser.add_argument("--read_length", type=int, default=150, help="Read length for ART")
    parser.add_argument("--insert_size", type=int, default=200, help="Mean insert size")
    parser.add_argument("--std_dev", type=int, default=10, help="Standard deviation of insert size")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    paternal_genome = args.paternal_fasta
    maternal_genome = args.maternal_fasta

    print(f"Simulating normal bulk reads at {args.normal_coverage}x coverage...")
    pat_normal_fa = os.path.join(args.output_dir, "paternal_normal.fa")
    mat_normal_fa = os.path.join(args.output_dir, "maternal_normal.fa")
    
    apply_cnvs_to_fasta(paternal_genome, pat_normal_fa, [args.het_cnp_bed], "pat")
    apply_cnvs_to_fasta(maternal_genome, mat_normal_fa, [args.het_cnp_bed], "mat")

    norm_pat_prefix = os.path.join(args.output_dir, "norm_pat_")
    norm_mat_prefix = os.path.join(args.output_dir, "norm_mat_")
    run_art(pat_normal_fa, args.normal_coverage / 2.0, args.read_length, args.insert_size, args.std_dev, norm_pat_prefix)
    run_art(mat_normal_fa, args.normal_coverage / 2.0, args.read_length, args.insert_size, args.std_dev, norm_mat_prefix)
    
    merge_fastqs(args.output_dir, "normal_bulk", [norm_pat_prefix, norm_mat_prefix])
    
    for f in [pat_normal_fa, mat_normal_fa]:
        if os.path.exists(f): os.remove(f)
    for p in [norm_pat_prefix, norm_mat_prefix]:
        for ext in ["1.fq", "2.fq"]:
            if os.path.exists(f"{p}{ext}"): os.remove(f"{p}{ext}")

    cell_vcfs = glob.glob(os.path.join(args.cell_dir, "*_snvs.vcf.gz"))
    if not cell_vcfs:
        cell_vcfs = glob.glob(os.path.join(args.cell_dir, "*", "*_snvs.vcf.gz"))

    num_cells = len(cell_vcfs)
    if num_cells == 0:
        raise ValueError("No *_snvs.vcf.gz files found in --cell_dir")

    print(f"Found {num_cells} cells. Simulating tumor reads and kindred clones...")
    
    # CALCULATE BOTH COVERAGES: Full for standalone clones, Fractional for the bulk mixture
    full_allele_cov = args.tumor_coverage / 2.0
    fraction_allele_cov = args.tumor_coverage / (2.0 * num_cells)
    
    bulk_part_prefixes_for_merge = []

    for vcf_file in cell_vcfs:
        cell_name = os.path.basename(vcf_file).replace("_snvs.vcf.gz", "")
        bed_file = vcf_file.replace("_snvs.vcf.gz", "_cnvs.bed")
        
        pat_snv_fa = os.path.join(args.output_dir, f"{cell_name}_pat_snv.fa")
        mat_snv_fa = os.path.join(args.output_dir, f"{cell_name}_mat_snv.fa")
        pat_final_fa = os.path.join(args.output_dir, f"{cell_name}_pat_final.fa")
        mat_final_fa = os.path.join(args.output_dir, f"{cell_name}_mat_final.fa")
        
        apply_bcftools_consensus(paternal_genome, vcf_file, pat_snv_fa, haplotype=1)
        apply_bcftools_consensus(maternal_genome, vcf_file, mat_snv_fa, haplotype=2)
        
        cnv_beds = [args.het_cnp_bed, bed_file]
        apply_cnvs_to_fasta(pat_snv_fa, pat_final_fa, cnv_beds, "pat")
        apply_cnvs_to_fasta(mat_snv_fa, mat_final_fa, cnv_beds, "mat")
        
        # ---------------------------------------------------------
        # 1. SIMULATE FULL DEPTH FOR STANDALONE CLONE FASTQ
        # ---------------------------------------------------------
        clone_pat_prefix = os.path.join(args.output_dir, f"{cell_name}_full_pat_")
        clone_mat_prefix = os.path.join(args.output_dir, f"{cell_name}_full_mat_")
        
        run_art(pat_final_fa, full_allele_cov, args.read_length, args.insert_size, args.std_dev, clone_pat_prefix)
        run_art(mat_final_fa, full_allele_cov, args.read_length, args.insert_size, args.std_dev, clone_mat_prefix)
        
        clone_out_prefix = f"{cell_name}_clone"
        merge_fastqs(args.output_dir, clone_out_prefix, [clone_pat_prefix, clone_mat_prefix])
        
        # ---------------------------------------------------------
        # 2. SIMULATE FRACTIONAL DEPTH FOR TUMOR BULK MIXTURE
        # ---------------------------------------------------------
        bulk_pat_prefix = os.path.join(args.output_dir, f"{cell_name}_frac_pat_")
        bulk_mat_prefix = os.path.join(args.output_dir, f"{cell_name}_frac_mat_")
        
        run_art(pat_final_fa, fraction_allele_cov, args.read_length, args.insert_size, args.std_dev, bulk_pat_prefix)
        run_art(mat_final_fa, fraction_allele_cov, args.read_length, args.insert_size, args.std_dev, bulk_mat_prefix)
        
        bulk_part_out_prefix = f"{cell_name}_bulk_part"
        merge_fastqs(args.output_dir, bulk_part_out_prefix, [bulk_pat_prefix, bulk_mat_prefix])
        
        bulk_part_prefixes_for_merge.append(os.path.join(args.output_dir, f"{bulk_part_out_prefix}_R"))
        
        # Cleanup temp cell FASTA files and intermediate split fastqs
        for f in [pat_snv_fa, mat_snv_fa, pat_final_fa, mat_final_fa]:
            if os.path.exists(f): os.remove(f)
        for p in [clone_pat_prefix, clone_mat_prefix, bulk_pat_prefix, bulk_mat_prefix]:
            for ext in ["1.fq", "2.fq"]:
                if os.path.exists(f"{p}{ext}"): os.remove(f"{p}{ext}")

    print("Merging generated fractional clone fastq files into match tumor bulk output...")
    r1_bulk = os.path.join(args.output_dir, "match_tumor_bulk_R1.fq")
    r2_bulk = os.path.join(args.output_dir, "match_tumor_bulk_R2.fq")
    
    with open(r1_bulk, "wb") as out1, open(r2_bulk, "wb") as out2:
        for prefix in bulk_part_prefixes_for_merge:
            r1_part = f"{prefix}1.fq"
            r2_part = f"{prefix}2.fq"
            if os.path.exists(r1_part):
                with open(r1_part, "rb") as f: out1.write(f.read())
                os.remove(r1_part) # Cleanup bulk parts after merging
            if os.path.exists(r2_part):
                with open(r2_part, "rb") as f: out2.write(f.read())
                os.remove(r2_part) # Cleanup bulk parts after merging

    print("Done! Match Tumor Bulk, Normal Bulk, and Full-Depth Clone FASTQs saved to:", args.output_dir)

if __name__ == "__main__":
    main()
