#!/bin/bash

set -e

# === INPUT ===
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --true-vcf) TRUE_VCF="$2"; shift ;;
        --test-dir) TEST_DIR="$2"; shift ;;
        --output) OUT_DIR="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

mkdir -p "$OUT_DIR"
echo "✅ True SNVs VCF: $TRUE_VCF"
echo "📁 Directory with test VCFs: $TEST_DIR"

# === Ensure true VCF is compressed and indexed ===
if [[ ! "$TRUE_VCF" =~ \.gz$ ]]; then
    echo "⚠️ Compressing $TRUE_VCF"
    bgzip -c "$TRUE_VCF" > "$TRUE_VCF.gz"
    TRUE_VCF="$TRUE_VCF.gz"
fi
# Add the -f flag to overwrite index if it already exists
[[ ! -f "$TRUE_VCF.tbi" ]] && bcftools index -f "$TRUE_VCF"

# === Loop through all test VCFs in the directory matching the wildcard ===
for TEST_VCF_DIR in $(find "$TEST_DIR" -type d); do
    for TEST_VCF in "$TEST_VCF_DIR"/*.vcf*; do
        TOOL_NAME=$(basename "$TEST_VCF" | cut -d. -f1)
        TOOL_OUT="$OUT_DIR/$TOOL_NAME"
        mkdir -p "$TOOL_OUT"

        echo ""
        echo "🔍 Benchmarking: $TOOL_NAME"

        # Compress/index if needed
        if [[ ! "$TEST_VCF" =~ \.gz$ ]]; then
            bgzip -c "$TEST_VCF" > "$TOOL_OUT/test.vcf.gz"
            TEST_VCF="$TOOL_OUT/test.vcf.gz"
        else
            cp "$TEST_VCF" "$TOOL_OUT/test.vcf.gz"
        fi
        [[ ! -f "$TEST_VCF.tbi" ]] && bcftools index -f "$TEST_VCF"

        # Intersect: True Positives
        bcftools isec -n=2 -w1 -Oz -o "$TOOL_OUT/true_positives.vcf.gz" "$TRUE_VCF" "$TEST_VCF"
        bcftools index "$TOOL_OUT/true_positives.vcf.gz"

        # Count
        TP=$(bcftools view -H "$TOOL_OUT/true_positives.vcf.gz" | wc -l)
        FN=$(bcftools view -H "$TRUE_VCF" | wc -l)
        FP=$(bcftools view -H "$TEST_VCF" | wc -l)
        FN=$((FN - TP))
        FP=$((FP - TP))

        # Metrics
        PREC=$(awk "BEGIN { if ($TP+$FP == 0) print 0; else print $TP / ($TP + $FP) }")
        RECALL=$(awk "BEGIN { if ($TP+$FN == 0) print 0; else print $TP / ($TP + $FN) }")
        F1=$(awk "BEGIN { if ($PREC+$RECALL == 0) print 0; else print 2 * $PREC * $RECALL / ($PREC + $RECALL) }")

        echo -e "$TOOL_NAME\tTP=$TP\tFP=$FP\tFN=$FN\tPrecision=$PREC\tRecall=$RECALL\tF1=$F1" | tee -a "$OUT_DIR/summary.tsv"
    done
done

echo ""
echo "🎉 Done! See summary in: $OUT_DIR/summary.tsv"
