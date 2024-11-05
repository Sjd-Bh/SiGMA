import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2

# Function to load SNV positions from multiple custom VCF files (maternal and paternal)
def load_combined_snv_positions(vcf_files):
    combined_snvs = set()
    for vcf_file in vcf_files:
        vcf_data = pd.read_csv(vcf_file, sep="\t", header=None, usecols=[0, 1], names=["chrom", "pos"])
        snv_positions = set(zip(vcf_data["chrom"], vcf_data["pos"]))
        combined_snvs.update(snv_positions)
    return combined_snvs

# Function to load SNV positions from ProSolo and SCcaller VCF files
def load_vcf_positions(vcf_file):
    vcf_data = pd.read_csv(vcf_file, sep="\t", comment='#', header=None, usecols=[0, 1], names=["chrom", "pos"])
    return set(zip(vcf_data["chrom"], vcf_data["pos"]))

# Function to calculate True Positives, False Positives, and False Negatives
def calculate_tp_fp_fn(detected_snvs_set, true_snvs_set):
    tp = len(detected_snvs_set & true_snvs_set)
    fp = len(detected_snvs_set - true_snvs_set)
    fn = len(true_snvs_set - detected_snvs_set)
    return tp, fp, fn

# Function to calculate Precision, Recall, and F1-Score
def calculate_metrics(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1_score

# Function to plot box plot of Precision and F1-Score
def plot_metrics_boxplot(data):
    plt.figure(figsize=(8, 5))
    sns.boxplot(x='Metric', y='Value', hue='Tool', data=data)
    plt.title("Precision and F1-Score for ProSolo and SCcaller")
    plt.ylabel("Score")
    plt.show()

# Function to plot Venn Diagram of detected SNVs by each tool
def plot_venn(prosolo_snvs_set, sccaller_snvs_set):
    plt.figure(figsize=(6, 6))
    venn2([prosolo_snvs_set, sccaller_snvs_set], ('ProSolo', 'SCcaller'))
    plt.title("Venn Diagram of Detected Variants")
    plt.show()

# Main function to load data, calculate metrics, and plot results
def main():
    # Load true SNV positions from the custom-format VCF file
    true_snv_files = ["maternal_snvs.vcf", "paternal_snvs.vcf"]  # Adjust paths as needed
    
    # Load combined true SNVs from both files
    true_snvs_set = load_combined_snv_positions(true_snv_files)

    # Load detected SNVs from ProSolo and SCcaller VCF files
    prosolo_snvs_set = load_vcf_positions("prosolo.vcf")        # Adjust path if needed
    sccaller_snvs_set = load_vcf_positions("sccaller.vcf")      # Adjust path if needed

    # Calculate TP, FP, and FN for ProSolo
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

# Run the main function
if __name__ == "__main__":
    main()
