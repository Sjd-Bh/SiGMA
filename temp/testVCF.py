# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 21:05:47 2024

@author: Sajedeh
"""

def main():
    # Load true SNV positions from the custom-format VCF file
    true_snv_files = ["../../test\\vc\\maternal_snvs_edit.vcf", 
                      "../../test\\vc\\paternal_snvs_edit.vcf"]  # Adjust paths as needed
    
    # Load combined true SNVs from both files
    true_snvs_set = load_combined_snv_positions(true_snv_files)
    gatk = load_vcf_positions("../../test\\vc\\bcf_calls.vcf") 
    # Load detected SNVs from ProSolo and SCcaller VCF files
    prosolo_snvs_set = load_vcf_positions("../../test\\vc\\prosolo_sim8_mda.vcf")        # Adjust path if needed
    sccaller_snvs_set = load_vcf_positions("../../test\\vc\\sccaller_sim2_mda.vcf")      # Adjust path if needed

    # Calculate TP, FP, and FN for ProSolo
    gatk_tp, gatk_fp, gatk_fn = calculate_tp_fp_fn(gatk, true_snvs_set)

    prosolo_tp, prosolo_fp, prosolo_fn = calculate_tp_fp_fn(prosolo_snvs_set, true_snvs_set)

    # Calculate TP, FP, and FN for SCcaller
    sccaller_tp, sccaller_fp, sccaller_fn = calculate_tp_fp_fn(sccaller_snvs_set, true_snvs_set)

    # Calculate metrics for ProSolo
    prosolo_precision, prosolo_recall, prosolo_f1 = calculate_metrics(prosolo_tp, prosolo_fp, prosolo_fn)

    # Calculate metrics for SCcaller
    sccaller_precision, sccaller_recall, sccaller_f1 = calculate_metrics(sccaller_tp, sccaller_fp, sccaller_fn)

    # Prepare data for plotting
    data = pd.DataFrame({
        'Tool': ['ProSolo', 'ProSolo', 'SCcaller', 'SCcaller'],
        'Metric': ['Precision', 'F1-score', 'Precision', 'F1-score'],
        'Value': [prosolo_precision, prosolo_f1, sccaller_precision, sccaller_f1]
    })

    # Plot metrics boxplot
    plot_metrics_boxplot(data)

    # Plot Venn Diagram
    plot_venn(prosolo_snvs_set, sccaller_snvs_set)
