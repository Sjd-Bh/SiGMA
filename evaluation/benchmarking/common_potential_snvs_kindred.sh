#!/bin/bash

# === Default Values ===
SC_DIR=""
BULK_NORMAL=""
BULK_KINDRED=""
OUTPUT_DIR="results"

# === Parse Arguments ===
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sc-dir) SC_DIR="$2"; shift ;;
        --bulk-normal) BULK_NORMAL="$2"; shift ;;
        --bulk-kindred) BULK_KINDRED="$2"; shift ;;
        --output) OUTPUT_DIR="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# === Validate Inputs ===
if [[ -z "$SC_DIR" || -z "$BULK_NORMAL" || -z "$BULK_KINDRED" ]]; then
    echo "Usage: $0 --sc-dir <path> --bulk-normal <vcf.gz> --bulk-kindred <vcf.gz> [--output <dir>]"
    exit 1
fi

# === Output dirs ===
TMP_DIR="$OUTPUT_DIR/tmp"
FINAL_DIR="$OUTPUT_DIR/final"
mkdir -p "$TMP_DIR" "$FINAL_DIR"

echo "📁 Single-cell VCFs from: $SC_DIR"
echo "🧬 Bulk normal VCF: $BULK_NORMAL"
echo "🧬 Kindred bulk VCF: $BULK_KINDRED"
echo "📂 Output directory: $OUTPUT_DIR"

# === 1. Compress and index all VCFs in SC_DIR if not already ===
#echo "📦 Compressing & indexing single-cell VCFs if needed..."
#for f in "$SC_DIR"/*.vcf; do
#    [[ -f "$f.gz" ]] || (bgzip -c "$f" > "$f.gz" && tabix -p vcf "$f.gz")
#done

# === 2. Merge single-cell VCFs ===
echo "🔄 Merging single-cell VCFs..."
bcftools merge --force-samples "$SC_DIR"/*.vcf.gz -Oz -o "$TMP_DIR/sc_merged.vcf.gz"
bcftools index "$TMP_DIR/sc_merged.vcf.gz"

# === 3. Remove variants found in bulk normal ===
echo "🚫 Removing variants found in bulk normal..."
bcftools isec -C -w1 -Oz -o "$TMP_DIR/sc_unique.vcf.gz" "$TMP_DIR/sc_merged.vcf.gz" "$BULK_NORMAL"
bcftools index "$TMP_DIR/sc_unique.vcf.gz"

# === 4. Keep variants seen in ≥2 cells ===
echo "🧮 Filtering shared variants (N_SAMPLES > 1)..."
bcftools view -i 'N_SAMPLES > 1' "$TMP_DIR/sc_unique.vcf.gz" -Oz -o "$FINAL_DIR/sc_shared.vcf.gz"
bcftools index "$FINAL_DIR/sc_shared.vcf.gz"

# === 5. Compare shared SNVs with kindred bulk ===
echo "🔬 Comparing with kindred bulk sample..."
bcftools isec "$FINAL_DIR/sc_shared.vcf.gz" "$BULK_KINDRED" -p "$OUTPUT_DIR/sc_vs_kindred"

# === 6. Rename final result and ensure valid compression ===
TRUE_SNV="$FINAL_DIR/true_snvs.vcf.gz"
if [[ -f "$OUTPUT_DIR/sc_vs_kindred/0003.vcf" ]]; then
    echo "✅ Extracting true SNVs..."
    mv "$OUTPUT_DIR/sc_vs_kindred/0003.vcf" "$FINAL_DIR/true_snvs.vcf"
    bgzip -f "$FINAL_DIR/true_snvs.vcf"
    bcftools index "$TRUE_SNV"
else
    echo "⚠️ No overlapping SNVs found. Skipping true SNV generation."
fi

# === Done ===
echo "✅ Benchmarking completed for single-cell samples."
echo "📄 Results in: $OUTPUT_DIR/"
