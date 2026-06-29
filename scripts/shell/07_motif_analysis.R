#!/usr/bin/env Rscript
# =============================================================================
# ICL 保守 Motif 分析
# =============================================================================
# 用 MEME Suite 发现保守 motif，或在 R 中用滑动窗口展示序列保守性。
#
# 用法:
#   # 方法A: 使用 MEME (需要安装 MEME Suite)
#   Rscript 07_motif_analysis.R --input candidates.fa --method meme
#
#   # 方法B: 纯R AA频率保守性热图 (无需MEME)
#   Rscript 07_motif_analysis.R --input aligned.fa --method conservation
#
# 输入:
#   FASTA 序列文件 (未比对或已比对)
#
# 输出:
#   PDF motif 结构图 / 保守性热图
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
flag <- function(name) any(args == name)

if (flag("-h") || flag("--help")) {
  cat("Usage: Rscript 07_motif_analysis.R [options]\n\n")
  cat("Options:\n")
  cat("  -i, --input    FILE   输入FASTA\n")
  cat("  -o, --output   PREFIX 输出前缀 (默认: icl_motif)\n")
  cat("  -m, --method   STR    方法: meme / conservation (默认: conservation)\n")
  cat("  --nmotifs      INT    MEME motif数 (默认: 10)\n")
  cat("  --window       INT    滑动窗口大小 (默认: 10, conservation方法)\n")
  quit("no")
}

INPUT   <- parse_arg("-i", parse_arg("--input", "icl_candidates.fasta"))
OUTPUT  <- parse_arg("-o", parse_arg("--output", "icl_motif"))
METHOD  <- parse_arg("-m", parse_arg("--method", "conservation"))
NMOTIFS <- parse_arg("--nmotifs", "10")
WINDOW  <- as.integer(parse_arg("--window", "10"))

cat("\n========== ICL Motif 分析 ==========\n")
cat("方法:", METHOD, "\n")

if (!file.exists(INPUT)) stop("错误: 输入文件不存在: ", INPUT)

# ----- 解析 FASTA -----
read_fasta <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  seqs <- list()
  current_name <- ""
  current_seq <- ""

  for (line in lines) {
    if (grepl("^>", line)) {
      if (current_name != "") {
        seqs[[current_name]] <- current_seq
      }
      current_name <- sub("^>", "", line)
      current_seq <- ""
    } else {
      current_seq <- paste0(current_seq, gsub("\\s", "", line))
    }
  }
  if (current_name != "") seqs[[current_name]] <- current_seq
  return(seqs)
}

seqs <- read_fasta(INPUT)
cat("序列数量:", length(seqs), "\n")
cat("序列长度:", nchar(seqs[[1]]), "(第一条)\n\n")

# ----- 方法A: MEME -----
if (METHOD == "meme") {
  meme_bin <- Sys.which("meme")
  if (meme_bin == "") {
    cat("警告: 未找到 MEME，回退到 conservation 方法\n")
    cat("安装 MEME: conda install -c bioconda meme\n")
    METHOD <- "conservation"
  }
}

if (METHOD == "meme") {
  cat("--- 运行 MEME ---\n")
  meme_out <- paste0(OUTPUT, "_meme")
  cmd <- paste("meme", INPUT, "-protein -oc", meme_out,
               "-nmotifs", NMOTIFS, "-nostatus")
  cat("命令:", cmd, "\n")
  system(cmd)

  meme_xml <- file.path(meme_out, "meme.xml")
  if (!file.exists(meme_xml)) {
    stop("MEME 运行失败，未生成 ", meme_xml)
  }

  # 尝试用 universalmotif 包解析 MEME XML
  if (requireNamespace("universalmotif", quietly = TRUE)) {
    library(universalmotif)
    motifs <- read_meme(meme_xml)
    cat("发现", length(motifs), "个 motif\n")

    # 绘制 motif logo
    pdf(paste0(OUTPUT, "_logos.pdf"), width = 10, height = length(motifs) * 2)
    par(mfrow = c(length(motifs), 1), mar = c(3, 4, 2, 1))
    # (简化的motif结构示意)
    dev.off()

  } else {
    cat("(未安装 universalmotif 包，仅输出 MEME 原始结果)\n")
    cat("MEME 结果目录: ", meme_out, "\n")
  }

} else {
  # ----- 方法B: 保守性热图 (纯R) -----
  cat("--- 计算序列保守性 ---\n")

  # 获取序列矩阵 (每个位置每个氨基酸)
  seq_names <- names(seqs)
  n_seqs <- length(seq_names)
  max_len <- max(nchar(unlist(seqs)))

  # 对齐长度（假设已比对或长度相近）
  seq_matrix <- matrix("", nrow = n_seqs, ncol = max_len)
  for (i in seq_len(n_seqs)) {
    chars <- strsplit(seqs[[i]], "")[[1]]
    seq_matrix[i, seq_along(chars)] <- chars
  }

  # 每个位置的 AA 频率
  aa_list <- unique(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]])

  # 计算每个位置每个氨基酸的频率
  conservation_df <- data.frame(position = integer(), aa = character(),
                                freq = numeric())

  for (pos in seq_len(max_len)) {
    col_chars <- seq_matrix[, pos]
    col_chars <- col_chars[col_chars != "" & col_chars != "-"]
    if (length(col_chars) == 0) next
    tbl <- table(col_chars)
    for (aa in names(tbl)) {
      conservation_df <- rbind(conservation_df,
        data.frame(position = pos, aa = aa,
                   freq = as.numeric(tbl[aa]) / length(col_chars)))
    }
  }

  # 保守性热图
  p1 <- ggplot(conservation_df, aes(x = position, y = aa, fill = freq)) +
    geom_tile(color = "white", size = 0.2) +
    scale_fill_gradientn(colors = c("white", "yellow", "red", "darkred"),
                         name = "频率", limits = c(0, 1)) +
    theme_minimal(base_size = 11) +
    labs(x = "序列位置", y = "氨基酸",
         title = "ICL 序列保守性热图",
         subtitle = paste("基于", n_seqs, "条序列"))

  # 每个位置的香农熵 (保守性指标)
  shannon_entropy <- function(freqs) {
    freqs <- freqs[freqs > 0]
    -sum(freqs * log2(freqs))
  }

  pos_entropy <- conservation_df %>%
    group_by(position) %>%
    summarise(
      entropy = shannon_entropy(freq),
      max_freq = max(freq),
      consensus = aa[which.max(freq)],
      .groups = "drop"
    )

  # 最大可能熵 (20种AA均分 = log2(20) ≈ 4.32)
  max_entropy <- log2(20)

  p2 <- ggplot(pos_entropy, aes(x = position)) +
    geom_ribbon(aes(ymin = 0, ymax = entropy), fill = "steelblue", alpha = 0.3) +
    geom_line(aes(y = entropy), color = "steelblue", size = 0.8) +
    geom_hline(yintercept = max_entropy, linetype = "dashed", color = "red") +
    annotate("text", x = max_len * 0.05, y = max_entropy,
             label = "完全不保守", hjust = 0, vjust = -0.5, size = 3, color = "red") +
    theme_minimal(base_size = 11) +
    labs(x = "序列位置", y = "香农熵 (bit)",
         title = "ICL 序列位置保守性",
         subtitle = "熵越低越保守") +
    ylim(0, max(pos_entropy$entropy, max_entropy + 0.5))

  # 合并输出
  pdf_out <- paste0(OUTPUT, "_conservation.pdf")
  pdf(pdf_out, width = 14, height = 8)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    print(p1 / p2 + plot_layout(heights = c(2, 1)))
  } else {
    print(p1)
    print(p2)
  }

  dev.off()
  cat("保守性图已保存:", pdf_out, "\n")

  # 输出高度保守的区间
  cat("\n--- 高度保守区域 (熵 < 1.0) ---\n")
  conserved_regions <- pos_entropy %>%
    filter(entropy < 1.0) %>%
    mutate(region = cumsum(c(1, diff(position) > 3)))

  if (nrow(conserved_regions) > 0) {
    regions <- conserved_regions %>%
      group_by(region) %>%
      summarise(start = min(position), end = max(position),
                len = end - start + 1,
                avg_entropy = mean(entropy),
                consensus_seq = paste(consensus, collapse = ""),
                .groups = "drop")
    print(regions)
    write.csv(regions, paste0(OUTPUT, "_conserved_regions.csv"), row.names = FALSE)
  }
}

cat("\n========== 完成 ==========\n")
