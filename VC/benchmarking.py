import argparse
import csv

def read_snvs_from_vcf(file_path):
    snvs = set()
    with open(file_path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue  # Skip headers and blank lines
            cols = line.strip().split()
            if len(cols) < 5:
                continue  # Skip malformed lines
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            snvs.add((chrom, pos, ref, alt))
    return snvs

def compute_metrics(truth_snvs, predicted_snvs):
    tp = len(truth_snvs & predicted_snvs)
    fp = len(predicted_snvs - truth_snvs)
    fn = len(truth_snvs - predicted_snvs)

    precision = tp / (tp + fp) if tp + fp > 0 else 0.0
    recall    = tp / (tp + fn) if tp + fn > 0 else 0.0
    f1        = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0.0

    return {
        "TruePositives": tp,
        "FalsePositives": fp,
        "FalseNegatives": fn,
        "Precision": round(precision, 4),
        "Recall": round(recall, 4),
        "F1Score": round(f1, 4)
    }

def write_metrics_to_csv(metrics, output_csv="benchmark_results.csv"):
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=metrics.keys())
        writer.writeheader()
        writer.writerow(metrics)

def main():
    parser = argparse.ArgumentParser(description="Benchmark SNV predictions against ground truth.")
    parser.add_argument("--ground-truth", required=True, help="Path to ground truth VCF")
    parser.add_argument("--pred", required=True, help="Path to predicted SNV VCF")
    parser.add_argument("--output", default="benchmark_results.csv", help="Output CSV file")
    args = parser.parse_args()

    truth_snvs = read_snvs_from_vcf(args.ground_truth)
    predicted_snvs = read_snvs_from_vcf(args.pred)

    metrics = compute_metrics(truth_snvs, predicted_snvs)
    write_metrics_to_csv(metrics, args.output)

    print("Benchmark results:")
    for k, v in metrics.items():
        print(f"{k}: {v}")

if __name__ == "__main__":
    main()

