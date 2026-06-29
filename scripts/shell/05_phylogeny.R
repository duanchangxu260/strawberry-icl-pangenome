#!/usr/bin/env Rscript
# =============================================================================
# ICL 基因家族系统发育树构建与可视化
# =============================================================================
# 输入: ICL 候选蛋白序列 FASTA (如 final_ICL_ogg.fa)
# 输出:
#   - 多序列比对 (aligned.fasta)
#   - 最大似然树 (icl_tree.nwk)
#   - 系统发育树图 PDF (icl_tree.pdf)
#
# 依赖: mafft/muscle, iqtree/fasttree (命令行)
#        R包: ggtree, treeio, ape, ggplot2, ggmsa (可选)
#
# 用法:
#   Rscript 05_phylogeny.R --input final_ICL_ogg.fa --output icl_tree
#   Rscript 05_phylogeny.R -i candidates.fa -o icl_tree --aligner muscle
#   Rscript 05_phylogeny.R -i candidates.fa --skip-alignment (用已有比对)
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(cowplot)
})

# ----- 参数解析 -----
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) > 0 && idx < length(args)) return(args[idx + 1])
  # handle --name=value format
  for (a in args) {
    if (grepl(paste0("^", name, "="), a)) return(sub(paste0("^", name, "="), "", a))
  }
  return(default)
}

flag <- function(name) any(args == name)

if (flag("-h") || flag("--help")) {
  cat("Usage: Rscript 05_phylogeny.R [options]\n\n")
  cat("Options:\n")
  cat("  -i, --input FILE    输入FASTA (默认: final_ICL_ogg.fa)\n")
  cat("  -o, --output PREFIX 输出前缀 (默认: icl_tree)\n")
  cat("  --aligner CMD       比对工具: mafft (默认) / muscle / clustalo\n")
  cat("  --tree-builder CMD  建树工具: iqtree (默认) / fasttree\n")
  cat("  --skip-alignment    使用已有比对 (输入为已比对FASTA)\n")
  cat("  --outgroup STR      外群标签 (默认自动选择拟南芥)\n")
  cat("  --threads INT       CPU线程数 (默认: 4)\n")
  quit("no")
}

INPUT    <- parse_arg("-i", parse_arg("--input", "final_ICL_ogg.fa"))
OUTPUT   <- parse_arg("-o", parse_arg("--output", "icl_tree"))
ALIGNER  <- parse_arg("--aligner", "mafft")
TREEBLD  <- parse_arg("--tree-builder", "iqtree")
OUTGROUP <- parse_arg("--outgroup", NULL)
THREADS  <- parse_arg("--threads", "4")
SKIP_ALN <- flag("--skip-alignment")

ALIGNED  <- paste0(OUTPUT, "_aligned.fasta")
TREEFILE <- paste0(OUTPUT, ".nwk")
PDFFILE  <- paste0(OUTPUT, ".pdf")

cat("\n========== ICL 系统发育树分析 ==========\n")
cat("输入序列:", INPUT, "\n")
cat("比对工具:", ALIGNER, "\n")
cat("建树工具:", TREEBLD, "\n")

if (!file.exists(INPUT)) {
  stop("错误: 输入文件不存在: ", INPUT)
}

# ----- Step 1: 多序列比对 -----
if (!SKIP_ALN) {
  cat("\n--- Step 1: 多序列比对 ---\n")

  aligner_cmd <- switch(ALIGNER,
    mafft   = paste("mafft --auto --thread", THREADS, INPUT, ">", ALIGNED),
    muscle  = paste("muscle -align", INPUT, "-output", ALIGNED),
    clustalo= paste("clustalo -i", INPUT, "-o", ALIGNED, "--threads", THREADS),
    stop("未知比对工具: ", ALIGNER)
  )

  cat("运行:", aligner_cmd, "\n")
  system(aligner_cmd)

  if (!file.exists(ALIGNED) || file.info(ALIGNED)$size == 0) {
    stop("比对失败，未生成输出文件")
  }
  cat("比对完成:", ALIGNED, "\n")
} else {
  ALIGNED <- INPUT
  cat("跳过比对，使用已有文件:", ALIGNED, "\n")
}

# ----- Step 2: 构建系统发育树 -----
cat("\n--- Step 2: 构建系统发育树 ---\n")

if (TREEBLD == "iqtree") {
  tree_cmd <- paste("iqtree -s", ALIGNED, "-nt", THREADS, "-pre", OUTPUT, "-quiet")
  cat("运行:", tree_cmd, "\n")
  system(tree_cmd)
  # IQ-TREE 输出 .treefile
  if (file.exists(paste0(OUTPUT, ".treefile"))) {
    file.copy(paste0(OUTPUT, ".treefile"), TREEFILE, overwrite = TRUE)
  }
} else if (TREEBLD == "fasttree") {
  tree_cmd <- paste("fasttree", ALIGNED, ">", TREEFILE)
  cat("运行:", tree_cmd, "\n")
  system(tree_cmd)
} else {
  stop("未知建树工具: ", TREEBLD)
}

if (!file.exists(TREEFILE) || file.info(TREEFILE)$size == 0) {
  stop("建树失败")
}
cat("建树完成:", TREEFILE, "\n")

# ----- Step 3: 读取并绘制系统发育树 -----
cat("\n--- Step 3: 绘制系统发育树 ---\n")

# 检测是否安装了 ggtree
has_ggtree <- requireNamespace("ggtree", quietly = TRUE)

if (has_ggtree && requireNamespace("treeio", quietly = TRUE)) {
  # ----- 方法A: ggtree (专业版) -----
  library(ggtree)
  library(treeio)

  tree <- read.newick(TREEFILE)

  # 自动检测外群 (拟南芥)
  if (is.null(OUTGROUP)) {
    tip_labels <- tree$tip.label
    ath_tips <- grep("[Aa]th|AT[0-9]|Arabidopsis", tip_labels, value = TRUE)
    if (length(ath_tips) > 0) {
      OUTGROUP <- ath_tips[1]
    }
  }

  if (!is.null(OUTGROUP) && OUTGROUP %in% tree$tip.label) {
    tree <- root(tree, outgroup = OUTGROUP, resolve.root = TRUE)
    cat("外群:", OUTGROUP, "\n")
  }

  # 主树
  p <- ggtree(tree, ladderize = TRUE, size = 0.8) +
    geom_tiplab(size = 3, offset = 0.02, align = FALSE) +
    geom_tippoint(size = 2, color = "steelblue") +
    theme_tree2() +
    labs(title = "ICL 基因家族系统发育树",
         subtitle = paste("Builder:", TREEBLD, "| Aligner:", ALIGNER))

  # 自举值 (如果有)
  if (!is.null(tree$node.label) && any(tree$node.label != "")) {
    p <- p + geom_nodelab(aes(label = ifelse(as.numeric(label) >= 70, label, "")),
                          size = 2.5, vjust = -0.5)
  }

  ggsave(PDFFILE, p, width = 10, height = max(6, length(tree$tip.label) * 0.4))
  cat("树图已保存:", PDFFILE, "\n")

} else {
  # ----- 方法B: ape 基础绘图 (无需 ggtree) -----
  cat("(未安装 ggtree，使用 ape 基础绘图)\n")

  tree <- read.tree(TREEFILE)

  if (!is.null(OUTGROUP) && OUTGROUP %in% tree$tip.label) {
    tree <- root(tree, outgroup = OUTGROUP, resolve.root = TRUE)
  }

  pdf(PDFFILE, width = 10, height = max(6, length(tree$tip.label) * 0.4))
  plot(tree, cex = 0.8, no.margin = TRUE)
  title(main = "ICL 基因家族系统发育树")
  add.scale.bar()
  # 添加自举值
  if (!is.null(tree$node.label) && any(tree$node.label != "")) {
    bs <- suppressWarnings(as.numeric(tree$node.label))
    bs[is.na(bs)] <- 0
    nodelabels(text = ifelse(bs >= 70, round(bs), ""),
               frame = "none", cex = 0.5, adj = c(1.2, -0.5))
  }
  dev.off()
  cat("树图已保存:", PDFFILE, "\n")
}

# ----- Step 4 (可选): MSA 可视化 -----
if (has_ggtree && requireNamespace("ggmsa", quietly = TRUE)) {
  cat("\n--- Step 4: MSA 可视化 ---\n")
  library(ggmsa)
  msa_pdf <- paste0(OUTPUT, "_msa.pdf")
  p_msa <- ggmsa(ALIGNED, start = 1, end = 200,
                 color = "Clustal", font = "DroidSansMono",
                 char_width = 0.5, seq_name = TRUE)
  ggsave(msa_pdf, p_msa, width = 14, height = max(4, length(tree$tip.label) * 0.3))
  cat("MSA 图已保存:", msa_pdf, "\n")
}

cat("\n========== 完成 ==========\n")
cat("输出文件:\n")
cat("  比对:", ALIGNED, "\n")
cat("  树文件:", TREEFILE, "\n")
cat("  树图:", PDFFILE, "\n")
