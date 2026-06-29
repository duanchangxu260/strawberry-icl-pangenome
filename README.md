# ICL Pan-Genome Analysis

**草莓（*Fragaria*）泛基因组异柠檬酸裂解酶（ICL）基因家族鉴定与分析流程。**

本仓库包含两套工具：
- **ICL 筛选** — HMMER → E-value 过滤 → 序列提取 → 拟南芥 BLAST → 直系同源组筛选
- **CDS 准备** — 泛基因组各物种 CDS 检测、批量提取、合并加前缀

---

## 依赖

```bash
# 必需
sudo apt install ncbi-blast+

# 推荐（可省略，有 awk 回退方案）
conda install -c bioconda seqkit gffread
```

不需要 Python。

---

## 目录结构

```
icl_pangenome_analysis/
├── README.md
├── LICENSE
└── scripts/
    ├── check_cds.sh              # 检测哪些物种缺 CDS
    ├── extract_cds_batch.sh      # 批量提取 CDS (gffread)
    ├── merge_rename_cds.sh       # 合并 CDS 并加物种前缀
    ├── 01_filter_hmm.sh          # HMMER E-value 筛选
    ├── 02_extract_sequences.sh   # 蛋白序列提取
    ├── 03_blast_arabidopsis.sh   # 拟南芥 BLAST + 过滤
    ├── 04_ogg_filter.R           # 直系同源组筛选 (R)
    ├── extract_fasta_by_ids.awk  # AWK FASTA 提取辅助
    └── run_pipeline.sh           # 一键运行
```

---

## 工作流程

### 前置: CDS 准备

```bash
# 1. 检查哪些基因组缺 CDS
bash scripts/check_cds.sh /data1/duanchangxu/Pan-Genome

# 2. 批量提取（需要每个物种目录下有基因组.fa + 注释.gff）
bash scripts/extract_cds_batch.sh /data1/duanchangxu/Pan-Genome

# 3. 合并并加物种前缀
bash scripts/merge_rename_cds.sh /data1/duanchangxu/Pan-Genome all_cds_merged.fasta
```

### ICL 基因家族筛选

**输入要求：** 每个物种子目录下有：
- `*.ICL.txt` — HMMER `--tblout` 输出
- `*.pep.fasta` 或 `*.fa` — 蛋白组

```bash
# 一键运行（四步）
bash scripts/run_pipeline.sh /mnt/e/ogg arabidopsis_proteome.fasta

# 或逐步运行
bash scripts/01_filter_hmm.sh /mnt/e/ogg 1e-5
bash scripts/02_extract_sequences.sh /mnt/e/ogg
bash scripts/03_blast_arabidopsis.sh /mnt/e/ogg ath.fa 1e-5 4 40 100
Rscript scripts/04_ogg_filter.R orthofinder --orthogroups Orthogroups.tsv --reference AT3G21720 -c best_hits.txt
```

### 步骤说明

| 步骤 | 脚本 | 输入 | 输出 |
|------|------|------|------|
| E-value 筛选 | `01_filter_hmm.sh` | `*.ICL.txt` | `*_filtered_ids.txt` |
| 序列提取 | `02_extract_sequences.sh` | ID列表 + 蛋白FASTA | `*_ICL_candidates.fasta` |
| 拟南芥 BLAST | `03_blast_arabidopsis.sh` | 候选序列 + ath.fa | `best_hits.txt` |
| OGG 筛选 | `04_ogg_filter.R` | best_hits + Orthogroups.tsv | `ogg_filtered_ids.txt` |

### 关键参数

- HMMER E-value: `≤ 1e-5`（可在命令中调整）
- BLAST E-value: `≤ 1e-5`, identity `≥ 40%`, 长度 `≥ 100aa`
- 拟南芥 ICL 参考基因: `AT3G21720`

### OGG 筛选三种方式

1. **OrthoFinder**（推荐）— 提供 `Orthogroups.tsv`，R 脚本自动匹配
2. **eggNOG** — 运行 `emapper.py` 后用 `04_ogg_filter.R` 匹配目标 OGG
3. **RBH** — 双向 BLAST 取互为最佳命中

---

## 常见问题

**ID 不匹配**: HMMER 输出的 ID 必须能匹配 FASTA 头的第一个空格分隔词。如果不一致，用 `seqkit grep -r` 或 `--match-by-prefix`。

**GFF/FASTA seqid 不一致**: 先 `grep "^>" genome.fa | head` 和 `cut -f1 annotation.gff | sort -u` 对比，用 `sed` 修正 GFF 第一列。

**路径含空格/中文**: gffread 对特殊字符敏感，建议用 `sed 's/ /_/g'` 统一重命名目录。

---

## 引用

- HMMER: Eddy SR (2011) *PLoS Comput Biol*
- BLAST+: Camacho C et al. (2009) *BMC Bioinformatics*
- OrthoFinder: Emms DM & Kelly S (2019) *Genome Biology*
- gffread: Pertea G & Pertea M (2020) *F1000Research*
