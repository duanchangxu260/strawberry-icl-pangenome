#!/usr/bin/env Rscript
# =============================================================================
# ICL 基因共线性分析 (Synteny)
# =============================================================================
# 基于 BLAST 或 MCScanX 结果，绘制 ICL 基因座周围的共线性图。
#
# 用法:
#   # 方法A: 使用 MCScanX 输出
#   Rscript 08_synteny_analysis.R --mcscanx output.collinearity
#
#   # 方法B: 简单的基因座对比图 (只需 BLAST + GFF)
#   Rscript 08_synteny_analysis.R --blast best_hits.txt --gff-dir gffs/
#
# 输出: PDF 共线性图
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  return(default)
}

if (any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript 08_synteny_analysis.R [options]\n\n")
  cat("Options:\n")
  cat("  --gff-dir    DIR    GFF3 目录 (用于基因座可视化)\n")
  cat("  --blast      FILE   BLAST 结果 (-outfmt 6)\n")
  cat("  --mcscanx    FILE   MCScanX collinearity 输出\n")
  cat("  -o, --output FILE   输出PDF (默认: icl_synteny.pdf)\n")
  cat("  --ref        STR    参考物种 (默认: Ath)\n")
  cat("  --flank      INT    侧翼基因数 (默认: 5)\n")
  quit("no")
}

GFF_DIR    <- parse_arg("--gff-dir", ".")
BLAST_FILE <- parse_arg("--blast", NULL)
MCSCANX    <- parse_arg("--mcscanx", NULL)
OUTPUT     <- parse_arg("-o", parse_arg("--output", "icl_synteny.pdf"))
REF_SPEC   <- parse_arg("--ref", "Ath")
FLANK      <- as.integer(parse_arg("--flank", "5"))

cat("\n========== ICL 共线性分析 ==========\n")

# ----- 解析 GFF，获取基因坐标 -----
parse_gff_genes <- function(gff_dir) {
  gff_files <- list.files(gff_dir, pattern = "\\.gff3?$", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)
  if (length(gff_files) == 0) return(NULL)

  all_genes <- list()
  for (f in gff_files) {
    lines <- readLines(f, warn = FALSE)
    lines <- lines[!grepl("^#", lines)]
    if (length(lines) == 0) next

    dat <- read.table(text = lines, sep = "\t", stringsAsFactors = FALSE,
                      quote = "", comment.char = "")
    colnames(dat) <- c("seqid", "source", "type", "start", "end",
                        "score", "strand", "phase", "attributes")

    gene_rows <- subset(dat, type == "gene")
    for (i in seq_len(nrow(gene_rows))) {
      attr <- gene_rows$attributes[i]
      gene_id <- regmatches(attr, regexpr("ID=([^;]+)", attr))
      gene_id <- sub("ID=", "", gene_id)
      gene_name <- regmatches(attr, regexpr("Name=([^;]+)", attr))
      gene_name <- sub("Name=", "", gene_name)

      all_genes[[length(all_genes) + 1]] <- data.frame(
        file     = basename(f),
        species  = gsub("(\\.gff3?|_ICL.*)$", "", basename(f)),
        seqid    = gene_rows$seqid[i],
        start    = gene_rows$start[i],
        end      = gene_rows$end[i],
        strand   = gene_rows$strand[i],
        gene_id  = ifelse(length(gene_id) > 0, gene_id, NA),
        gene_name= ifelse(length(gene_name) > 0, gene_name, NA),
        stringsAsFactors = FALSE
      )
    }
  }
  return(do.call(rbind, all_genes))
}

genes <- parse_gff_genes(GFF_DIR)

# ----- 方法A: MCScanX 结果 -----
if (!is.null(MCSCANX) && file.exists(MCSCANX)) {
  cat("解析 MCScanX 结果:", MCSCANX, "\n")

  # MCScanX collinearity 格式: 每行两个基因
  col_data <- read.table(MCSCANX, stringsAsFactors = FALSE, fill = TRUE)
  cat("共线性基因对数:", nrow(col_data), "\n")

  # 简化可视化: 显示 ICL 基因周围的共线性区块
  # 这里做一个简单的共线性连线图

  # 筛选含 ICL 的区块
  icl_blocks <- col_data[grepl("ICL|icl|AT3G21720", apply(col_data, 1, paste, collapse=" ")), ]
  cat("含 ICL 的区块:", nrow(icl_blocks), "\n")

  pdf(OUTPUT, width = 10, height = 6)
  plot(1, type = "n", xlim = c(0, 10), ylim = c(0, 10),
       xlab = "", ylab = "", axes = FALSE,
       main = "ICL 共线性区块 (MCScanX)")

  if (nrow(icl_blocks) > 0) {
    n_blocks <- min(nrow(icl_blocks), 20)
    for (i in seq_len(n_blocks)) {
      y <- 9 - i * (8 / n_blocks)
      lines(c(1, 9), c(y, y), lwd = 2, col = scales::hue_pal()(n_blocks)[i])
      text(0.5, y, icl_blocks[i, 1], cex = 0.6, adj = 1)
      text(9.5, y, icl_blocks[i, 2], cex = 0.6, adj = 0)
    }
  } else {
    text(5, 5, "未找到含 ICL 的共线性区块\n请检查 MCScanX 输入数据", cex = 1.2)
  }
  dev.off()

# ----- 方法B: 基因座对比图 -----
} else if (!is.null(genes) && nrow(genes) > 0) {
  cat("使用基因坐标绘制基因座图\n")

  # 找出 ICL 基因
  icl_genes <- genes[grepl("ICL|icl|AT3G21720", genes$gene_name) |
                     grepl("ICL|icl|AT3G21720", genes$gene_id), ]

  if (nrow(icl_genes) == 0) {
    # 如果没找到明确标注的 ICL，直接用所有基因
    icl_genes <- genes
    cat("注意: 未找到明确标注为 ICL 的基因，显示所有基因\n")
  } else {
    cat("找到", nrow(icl_genes), "个 ICL 基因\n")
  }

  pdf(OUTPUT, width = 12, height = max(4, nrow(icl_genes) * 1.2))

  # 按物种排列
  icl_genes <- icl_genes[order(icl_genes$species), ]
  icl_genes$y <- seq_len(nrow(icl_genes))

  p <- ggplot(icl_genes) +
    # 染色体/骨架
    geom_segment(aes(x = start - 5000, xend = end + 5000,
                     y = y, yend = y),
                 linewidth = 0.5, color = "grey70") +
    # 基因箭头
    geom_segment(aes(x = start, xend = end,
                     y = y, yend = y,
                     color = strand),
                 linewidth = 4, arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    # 基因标签
    geom_text(aes(x = (start + end) / 2, y = y + 0.5,
                  label = paste(gene_name, gene_id, sep = "\n")),
              size = 2.5, vjust = 0, lineheight = 0.8) +
    scale_color_manual(values = c("+" = "steelblue", "-" = "coral"),
                       name = "链") +
    scale_y_continuous(breaks = icl_genes$y,
                       labels = paste(icl_genes$species, icl_genes$seqid, sep = "\n")) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 7, lineheight = 0.8),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "基因组坐标 (bp)", y = "",
         title = "ICL 基因座对比图",
         subtitle = paste("显示", nrow(icl_genes), "个基因 | 侧翼 ±5kb"))

  print(p)
  dev.off()

} else {
  cat("错误: 既没有 MCScanX 结果也没有 GFF 数据\n")
  cat("请提供 --gff-dir 或 --mcscanx\n")
  quit(status = 1)
}

cat("共线性图已保存:", OUTPUT, "\n")
cat("\n========== 完成 ==========\n")
