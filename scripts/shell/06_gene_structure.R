#!/usr/bin/env Rscript
# =============================================================================
# ICL 基因结构图（外显子-内含子）
# =============================================================================
# 从 GFF3 文件提取 ICL 基因的 exon/intron 结构并绘制对比图。
#
# 用法:
#   Rscript 06_gene_structure.R --gff-dir test_data --gene-ids genes.txt
#   Rscript 06_gene_structure.R -d gff_files/ -g icl_genes.txt -o icl_structure.pdf
#
# 输入:
#   GFF3 文件（每个物种一个，放在 --gff-dir 或当前目录）
#   基因 ID 列表（每行一个 gene_id，对应 GFF 中 mRNA 的 Parent）
#
# 输出:
#   基因结构对比 PDF
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

if (any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript 06_gene_structure.R [options]\n\n")
  cat("Options:\n")
  cat("  -d, --gff-dir  DIR   GFF3 文件目录 (默认: .)\n")
  cat("  -g, --gene-ids FILE  基因 ID 列表 (每行一个)\n")
  cat("  -o, --output    FILE  输出PDF (默认: icl_gene_structure.pdf)\n")
  cat("  --width         NUM   PDF宽度 (默认: 10)\n")
  cat("  --height        NUM   PDF高度 (默认: 6)\n")
  quit("no")
}

GFF_DIR  <- parse_arg("-d", parse_arg("--gff-dir", "."))
GENE_IDS <- parse_arg("-g", parse_arg("--gene-ids", NULL))
OUTPUT   <- parse_arg("-o", parse_arg("--output", "icl_gene_structure.pdf"))
FIG_W    <- as.numeric(parse_arg("--width", "10"))
FIG_H    <- as.numeric(parse_arg("--height", "6"))

cat("\n========== ICL 基因结构图 ==========\n")

# ----- Step 1: 查找所有 GFF3 文件 -----
gff_files <- list.files(GFF_DIR, pattern = "\\.gff3?$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)

if (length(gff_files) == 0) {
  stop("在 ", GFF_DIR, " 中未找到 GFF3 文件")
}
cat("找到", length(gff_files), "个 GFF3 文件\n")

# ----- Step 2: 解析 GFF3，提取 exon 坐标 -----
parse_gff_exons <- function(gff_path) {
  if (!file.exists(gff_path)) return(NULL)

  # 跳过注释行读取 GFF
  lines <- readLines(gff_path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]

  if (length(lines) == 0) return(NULL)

  dat <- read.table(text = lines, sep = "\t", stringsAsFactors = FALSE,
                    quote = "", comment.char = "")
  colnames(dat) <- c("seqid", "source", "type", "start", "end",
                      "score", "strand", "phase", "attributes")

  # 过滤 exon 行
  exons <- subset(dat, type == "exon")

  if (nrow(exons) == 0) return(NULL)

  # 提取基因名和转录本名
  extract_attr <- function(attr_str, key) {
    m <- regmatches(attr_str, regexpr(paste0(key, "=([^;]+)"), attr_str))
    if (length(m) == 0) return(NA)
    sub(paste0(key, "="), "", m)
  }

  exons$Parent  <- sapply(exons$attributes, extract_attr, "Parent")
  exons$gene_id <- sapply(exons$attributes, extract_attr, "ID")
  exons$gene_id <- ifelse(is.na(exons$gene_id),
                          sapply(exons$attributes, extract_attr, "gene_id"),
                          exons$gene_id)

  # 从文件名推断物种名
  exons$species <- gsub("(\\.gff3?|_ICL.*)$", "", basename(gff_path))
  exons$species <- gsub("[._]ICL$", "", exons$species)

  return(exons)
}

cat("解析 GFF 文件中...\n")
all_exons <- do.call(rbind, lapply(gff_files, function(f) {
  res <- parse_gff_exons(f)
  if (!is.null(res)) res$file <- basename(f)
  return(res)
}))

if (is.null(all_exons) || nrow(all_exons) == 0) {
  stop("未能从任何 GFF 中提取到 exon 信息")
}

cat("共提取到", nrow(all_exons), "个 exon\n")

# ----- Step 3: 标准化为相对坐标 -----
# 对每个基因，用第一个 exon 的 start 作为 0
all_exons <- all_exons %>%
  group_by(file, Parent) %>%
  mutate(
    gene_start = min(start),
    gene_end   = max(end),
    rel_start  = start - gene_start,
    rel_end    = end - gene_start,
    gene_len   = gene_end - gene_start
  ) %>%
  ungroup()

# ----- Step 4: 绘制基因结构图 -----
cat("绘图...\n")

# 尝试用 gggenes；如果没装就回退到 ggplot2 手动绘制
has_gggenes <- requireNamespace("gggenes", quietly = TRUE)

if (has_gggenes) {
  library(gggenes)

  p <- ggplot(all_exons, aes(xmin = rel_start, xmax = rel_end,
                              y = file, forward = strand == "+")) +
    geom_gene_arrow(fill = "steelblue", color = "grey30", size = 0.5,
                    arrowhead_height = unit(4, "mm"),
                    arrowhead_width  = unit(3, "mm")) +
    geom_subgene_arrow(aes(xsubmin = rel_start, xsubmax = rel_end,
                           fill = "CDS"), color = NA) +
    facet_wrap(~ file, scales = "free", ncol = 1) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "相对位置 (bp)", y = "",
         title = "ICL 基因结构 (Exon-Intron)",
         subtitle = paste("共", length(unique(all_exons$file)), "个物种"))

} else {
  # 回退方案: ggplot2 手动绘制
  cat("(未安装 gggenes，使用 ggplot2 手动绘制)\n")

  # 给每个基因赋予唯一的 y 位置
  gene_labels <- unique(all_exons$file)
  all_exons$y_pos <- match(all_exons$file, gene_labels)

  # 区分正负链
  all_exons$dir <- ifelse(all_exons$strand == "+", 1, -1)

  p <- ggplot(all_exons) +
    # 基因骨架线
    geom_segment(aes(x = 0, xend = gene_len, y = y_pos, yend = y_pos),
                 linewidth = 1.2, color = "grey40") +
    # exon 矩形
    geom_rect(aes(xmin = rel_start, xmax = rel_end,
                  ymin = y_pos - 0.3, ymax = y_pos + 0.3),
              fill = "steelblue", color = "grey20", size = 0.3) +
    scale_y_continuous(breaks = seq_along(gene_labels),
                       labels = gene_labels) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "相对位置 (bp)", y = "",
         title = "ICL 基因结构 (Exon-Intron)",
         subtitle = paste("蓝色矩形 = exon, 灰色线 = intron | 共",
                          length(gene_labels), "个基因"))
}

ggsave(OUTPUT, p, width = FIG_W, height = FIG_H)
cat("基因结构图已保存:", OUTPUT, "\n")

# ----- Step 5: 输出统计表 -----
cat("\n--- 基因结构统计 ---\n")
stats <- all_exons %>%
  group_by(file, Parent, species) %>%
  summarise(
    exon_count = n(),
    total_cds_len = sum(end - start + 1),
    gene_span = max(gene_len),
    strand = first(strand),
    .groups = "drop"
  )

print(as.data.frame(stats))

# 保存统计表
write.csv(stats, sub("\\.pdf$", "_stats.csv", OUTPUT), row.names = FALSE)
cat("\n统计表已保存:", sub("\\.pdf$", "_stats.csv", OUTPUT), "\n")
cat("\n========== 完成 ==========\n")
