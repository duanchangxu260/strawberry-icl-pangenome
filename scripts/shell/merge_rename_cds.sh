#!/bin/bash
# =============================================================================
# 合并所有物种 CDS 并添加物种前缀（为下游泛基因组分析准备）
# =============================================================================
# 将每个物种子目录下的 cds.fasta 合并，并用 seqkit 给每条序列加物种名前缀。
# 无 seqkit 时用 awk 替代。
#
# 用法:
#   bash merge_rename_cds.sh <pan_genome_dir> [output_file]
#   bash merge_rename_cds.sh /data1/duanchangxu/Pan-Genome all_cds.fasta
# =============================================================================

set -euo pipefail

WORK_DIR="${1:-.}"
OUTPUT="${2:-all_cds_merged.fasta}"

if [ ! -d "$WORK_DIR" ]; then
    echo "错误: 目录不存在: $WORK_DIR"
    exit 1
fi

echo "== 合并 CDS 序列并添加物种前缀 =="
echo ""

tmp_file=$(mktemp)
trap 'rm -f "$tmp_file"' EXIT

count=0
species_count=0

for genome_dir in "$WORK_DIR"/*/; do
    genome_name=$(basename "$genome_dir")

    # 查找 CDS 文件
    cds_file=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.cds.fa" -o -name "*.cds.fasta" -o -name "*.cds.fna" \
        -o -name "*.cds.ffn" -o -name "*.cds" -o -name "cds.fasta" \) \
        | head -1)

    if [ -z "$cds_file" ]; then
        echo "⊙ $genome_name — 无 CDS，跳过"
        continue
    fi

    # 用物种名做前缀（把空格和特殊字符替换为下划线）
    prefix=$(echo "$genome_name" | sed 's/[^a-zA-Z0-9]/_/g' | sed 's/__*/_/g' | sed 's/_$//')

    seqs=$(grep -c "^>" "$cds_file" 2>/dev/null || echo 0)

    if command -v seqkit &> /dev/null; then
        seqkit rename -p "^($|.*)" -r "${prefix}_"'$1' "$cds_file" >> "$tmp_file"
    else
        # awk 回退: 给 > 行加前缀
        awk -v pfx="${prefix}_" '/^>/ { sub(/^>/, ">" pfx) } { print }' "$cds_file" >> "$tmp_file"
    fi

    echo "✓ $genome_name — $seqs 条序列 (前缀: ${prefix}_)"
    count=$((count + seqs))
    species_count=$((species_count + 1))
done

mv "$tmp_file" "$OUTPUT"

echo ""
echo "== 完成 =="
echo "物种数: $species_count"
echo "总 CDS 序列: $count"
echo "输出文件: $OUTPUT"
