#!/bin/bash
# =============================================================================
# 检测泛基因组目录下哪些物种缺少 CDS 序列文件
# =============================================================================
# 用法:
#   bash check_cds.sh <pan_genome_dir>
#   bash check_cds.sh /data1/duanchangxu/Pan-Genome
# =============================================================================

set -euo pipefail

WORK_DIR="${1:-.}"

if [ ! -d "$WORK_DIR" ]; then
    echo "错误: 目录不存在: $WORK_DIR"
    exit 1
fi

echo "== 检查 CDS 文件状态 =="
echo ""

missing=0
present=0

for genome_dir in "$WORK_DIR"/*/; do
    genome_name=$(basename "$genome_dir")

    # 查找 CDS 文件（支持多种扩展名）
    cds_file=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.cds.fa" -o -name "*.cds.fasta" -o -name "*.cds.fna" \
        -o -name "*.cds.ffn" -o -name "*.cds" -o -name "cds.fasta" \) \
        | head -1)

    if [ -n "$cds_file" ]; then
        count=$(grep -c "^>" "$cds_file" 2>/dev/null || echo "?")
        echo "✓ $genome_name — $(basename "$cds_file") ($count 条序列)"
        present=$((present + 1))
    else
        # 检查是否有基因组+注释可以提取
        has_genome=$(find "$genome_dir" -maxdepth 1 -type f \
            \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
            ! -name "*cds*" | head -1)
        has_gff=$(find "$genome_dir" -maxdepth 1 -type f \
            \( -name "*.gff" -o -name "*.gff3" -o -name "*.gtf" \) | head -1)

        if [ -n "$has_genome" ] && [ -n "$has_gff" ]; then
            echo "✗ $genome_name — 无 CDS（可提取: 有基因组+注释）"
        else
            echo "✗ $genome_name — 无 CDS（且缺基因组或注释，无法提取）"
        fi
        missing=$((missing + 1))
    fi
done

echo ""
echo "== 汇总: $present 个有CDS, $missing 个缺失 =="
