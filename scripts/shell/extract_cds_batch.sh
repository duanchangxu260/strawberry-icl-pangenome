#!/bin/bash
# =============================================================================
# 批量从基因组 + GFF 注释中提取 CDS 序列（使用 gffread）
# =============================================================================
# 遍历泛基因组目录下每个物种子目录，若缺少 CDS 文件则自动用 gffread 提取。
#
# 前置条件:
#   每个物种目录下需要有:
#     - 基因组 FASTA (*.fa / *.fasta / *.fna)
#     - 注释文件 (*.gff / *.gff3 / *.gtf)
#
# 用法:
#   bash extract_cds_batch.sh <pan_genome_dir>
#   bash extract_cds_batch.sh /data1/duanchangxu/Pan-Genome
#
# 来源: DeepSeek ICL 泛基因组分析流程
# =============================================================================

set -euo pipefail

WORK_DIR="${1:-.}"

if [ ! -d "$WORK_DIR" ]; then
    echo "错误: 目录不存在: $WORK_DIR"
    exit 1
fi

if ! command -v gffread &> /dev/null; then
    echo "错误: 未安装 gffread"
    echo "安装: conda install -c bioconda gffread"
    echo " 或:  sudo apt install gffread"
    exit 1
fi

echo "== 批量提取 CDS 序列 =="
echo ""

extracted=0
skipped=0
failed=0

for genome_dir in "$WORK_DIR"/*/; do
    genome_name=$(basename "$genome_dir")

    # 检查是否已有 CDS
    existing=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.cds.fa" -o -name "*.cds.fasta" -o -name "*.cds.fna" \
        -o -name "*.cds.ffn" -o -name "*.cds" -o -name "cds.fasta" \) \
        | head -1)

    if [ -n "$existing" ]; then
        echo "⊙ $genome_name — 已有 CDS，跳过"
        skipped=$((skipped + 1))
        continue
    fi

    # 查找基因组和注释文件
    GENOME=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
        ! -name "*cds*" ! -name "*pep*" ! -name "*protein*" | head -1)
    GFF=$(find "$genome_dir" -maxdepth 1 -type f \
        \( -name "*.gff" -o -name "*.gff3" -o -name "*.gtf" \) | head -1)

    if [ -z "$GENOME" ] || [ -z "$GFF" ]; then
        echo "✗ $genome_name — 缺少基因组或注释文件，跳过"
        failed=$((failed + 1))
        continue
    fi

    # 检查 GFF 中是否有 CDS 特征
    cds_features=$(grep -ci "[^#].*CDS" "$GFF" 2>/dev/null || echo 0)
    if [ "$cds_features" -eq 0 ]; then
        echo "✗ $genome_name — GFF 中无 CDS 特征，尝试 -g 提取"
        # 如果基因序列文件本身就是 gene.fasta (cDNA), 直接用
        gene_fa=$(find "$genome_dir" -maxdepth 1 -type f -name "*.gene.fasta" | head -1)
        if [ -n "$gene_fa" ]; then
            echo "  → 发现 gene.fasta, 直接作为 CDS 使用"
            cp "$gene_fa" "${genome_dir}/cds.fasta"
            extracted=$((extracted + 1))
            continue
        fi
    fi

    # 用 gffread 提取 CDS
    OUTPUT="${genome_dir}/cds.fasta"
    echo "→ $genome_name 提取中..."

    if gffread -g "$GENOME" -x "$OUTPUT" "$GFF" 2>/dev/null; then
        count=$(grep -c "^>" "$OUTPUT" 2>/dev/null || echo 0)
        echo "  ✓ 提取完成: $count 条 CDS"
        extracted=$((extracted + 1))
    else
        # 如果 gffread 因 seqid 不匹配失败，尝试用 gene.fasta 代替
        gene_fa=$(find "$genome_dir" -maxdepth 1 -type f -name "*.gene.fasta" | head -1)
        if [ -n "$gene_fa" ] && [ ! -f "$OUTPUT" ]; then
            echo "  ⚠ gffread 失败，回退使用 gene.fasta 作为 CDS"
            cp "$gene_fa" "$OUTPUT"
            count=$(grep -c "^>" "$OUTPUT" 2>/dev/null || echo 0)
            echo "  ✓ 复制完成: $count 条 CDS"
            extracted=$((extracted + 1))
        else
            echo "  ✗ 提取失败"
            failed=$((failed + 1))
            rm -f "$OUTPUT"
        fi
    fi
done

echo ""
echo "== 完成: $extracted 个新提取, $skipped 个跳过, $failed 个失败 =="
echo ""
echo "提示: 提取的 CDS 文件名为 cds.fasta，建议用 merge_rename_cds.sh 合并"

if [ "$failed" -gt 0 ]; then
    echo ""
    echo "常见失败原因:"
    echo "  1. GFF 的 seqid 与 FASTA 头不匹配 → 用 check_gff_fasta_match.sh 排查"
    echo "  2. GFF 中缺少 mRNA/CDS 行 → 检查注释完整性"
    echo "  3. 路径含空格/特殊字符 → 目录名中的单引号 ' 可能导致问题"
fi
