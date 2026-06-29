#!/usr/bin/env Rscript
# =============================================================================
# Step 4: Orthogroup (OGG) filtering — R implementation
# =============================================================================
# Filters ICL candidate sequences by orthogroup membership.
#
# Method A: Use OrthoFinder Orthogroups.tsv
#   Finds the orthogroup containing the reference gene (e.g. Arabidopsis ICL),
#   extracts all genes in that OG, intersects with candidate list.
#
# Method B: Reciprocal Best Hit (RBH) from bidirectional BLAST
#   Retains query-target pairs that are reciprocal best BLAST hits.
#
# Usage:
#   # Method A: OrthoFinder
#   Rscript 04_ogg_filter.R orthofinder \
#       --orthogroups Orthogroups.tsv \
#       --reference AT3G21720 \
#       --candidates best_hits.txt \
#       --output ogg_keep.txt
#
#   # Method B: Reciprocal Best Hit
#   Rscript 04_ogg_filter.R rbh \
#       --forward best_hits.txt \
#       --reverse ath_vs_cand.txt \
#       --candidates best_hits.txt \
#       --output rbh_keep.txt
#
# Input:
#   OrthoFinder: Orthogroups.tsv, candidate ID list / best_hits.txt
#   RBH: Forward and reverse BLAST results (outfmt 6)
#
# Output:
#   Filtered sequence ID list (one per line)
# =============================================================================

library(methods)  # Explicitly load for Rscript compatibility

# ----- Parse command-line arguments -----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript 04_ogg_filter.R <method> [options]\n\n")
  cat("Methods:\n")
  cat("  orthofinder  Filter using OrthoFinder Orthogroups.tsv\n")
  cat("  rbh          Reciprocal Best Hit filtering\n\n")
  cat("OrthoFinder options:\n")
  cat("  --orthogroups FILE   Orthogroups.tsv path\n")
  cat("  --reference   GENE   Reference gene ID (e.g., AT3G21720)\n")
  cat("  --candidates  FILE   Candidate ID list or BLAST best_hits.txt\n")
  cat("  --output      FILE   Output file (default: ogg_keep.txt)\n\n")
  cat("RBH options:\n")
  cat("  --forward     FILE   Forward BLAST results (candidates vs ref)\n")
  cat("  --reverse     FILE   Reverse BLAST results (ref vs candidates)\n")
  cat("  --candidates  FILE   Candidate file for final intersection\n")
  cat("  --output      FILE   Output file (default: rbh_keep.txt)\n")
  quit(status = 1)
}

method <- args[1]

# Parse named arguments
parse_arg <- function(name, default = NULL) {
  idx <- which(args == name)
  if (length(idx) > 0 && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default)
}

# ----- Method A: OrthoFinder filtering -----
method_orthofinder <- function() {
  orthogroups_file <- parse_arg("--orthogroups")
  reference_gene   <- parse_arg("--reference")
  candidates_file  <- parse_arg("--candidates")
  output_file      <- parse_arg("--output", "ogg_keep.txt")

  if (is.null(orthogroups_file) || is.null(reference_gene)) {
    stop("--orthogroups and --reference are required for orthofinder method")
  }

  cat("=== Method: OrthoFinder Orthogroup Filtering ===\n")
  cat("Orthogroups file:", orthogroups_file, "\n")
  cat("Reference gene:", reference_gene, "\n")

  # Read Orthogroups.tsv
  ortho <- read.delim(orthogroups_file, header = TRUE, stringsAsFactors = FALSE,
                      check.names = FALSE)
  cat(sprintf("Loaded %d orthogroups\n", nrow(ortho)))

  # Find the orthogroup containing the reference gene
  target_og <- NULL
  for (i in seq_len(nrow(ortho))) {
    for (col in colnames(ortho)[-1]) {  # Skip Orthogroup column
      genes_str <- as.character(ortho[i, col])
      if (grepl(reference_gene, genes_str, fixed = TRUE)) {
        target_og <- ortho[i, "Orthogroup"]
        break
      }
    }
    if (!is.null(target_og)) break
  }

  if (is.null(target_og)) {
    stop(sprintf("Reference gene '%s' not found in Orthogroups.tsv", reference_gene))
  }

  cat(sprintf("Found reference gene '%s' in orthogroup: %s\n",
              reference_gene, target_og))

  # Extract all genes in this orthogroup
  og_row <- ortho[ortho$Orthogroup == target_og, , drop = FALSE]
  genes_in_og <- character()

  for (col in colnames(ortho)[-1]) {
    cell <- as.character(og_row[1, col])
    if (!is.na(cell) && nchar(cell) > 0) {
      genes <- trimws(unlist(strsplit(cell, ",")))
      genes_in_og <- c(genes_in_og, genes)
    }
  }
  genes_in_og <- unique(genes_in_og[genes_in_og != ""])

  cat(sprintf("Total genes in orthogroup: %d\n", length(genes_in_og)))

  # Read candidate IDs
  if (!is.null(candidates_file) && file.exists(candidates_file)) {
    # Detect if it's a BLAST best_hits.txt (tab-separated) or a plain ID list
    first_line <- readLines(candidates_file, n = 1)
    if (grepl("\t", first_line)) {
      candidates <- unique(read.table(candidates_file, stringsAsFactors = FALSE)[, 1])
    } else {
      candidates <- unique(readLines(candidates_file))
    }
    candidates <- candidates[candidates != ""]

    cat(sprintf("Candidate sequences: %d\n", length(candidates)))

    # Intersect
    keep <- intersect(candidates, genes_in_og)
    cat(sprintf("After intersection with candidates: %d\n", length(keep)))
  } else {
    keep <- genes_in_og
    cat("No candidates file provided, outputting all genes in OG.\n")
  }

  # Write output
  writeLines(keep, output_file)
  cat(sprintf("\nFiltered IDs written to: %s\n", output_file))
  cat(sprintf("Total: %d sequences\n", length(keep)))
}

# ----- Method B: Reciprocal Best Hit (RBH) -----
method_rbh <- function() {
  forward_file  <- parse_arg("--forward")
  reverse_file  <- parse_arg("--reverse")
  candidates_file <- parse_arg("--candidates")
  output_file   <- parse_arg("--output", "rbh_keep.txt")

  if (is.null(forward_file) || is.null(reverse_file)) {
    stop("--forward and --reverse are required for rbh method")
  }

  cat("=== Method: Reciprocal Best Hit (RBH) ===\n")
  cat("Forward BLAST:", forward_file, "\n")
  cat("Reverse BLAST:", reverse_file, "\n")

  # Parse forward BLAST: best hit per query (lowest E-value in col 11 or 4)
  parse_best <- function(filepath, eval_col = NULL) {
    dat <- read.table(filepath, stringsAsFactors = FALSE, fill = TRUE,
                      comment.char = "", sep = "\t")
    # Auto-detect E-value column
    if (is.null(eval_col)) {
      # Try common positions (0-indexed): col 10 (standard), col 4 (simplified)
      if (ncol(dat) >= 11 && is.numeric(dat[, 11])) {
        eval_col <- 11
      } else if (ncol(dat) >= 5 && is.numeric(dat[, 5])) {
        eval_col <- 5
      } else {
        stop("Cannot auto-detect E-value column in BLAST output")
      }
    }
    # For each query (col 1), keep row with min E-value
    best <- do.call(rbind, by(dat, dat[, 1], function(x) x[which.min(x[, eval_col]), ]))
    best <- best[, c(1, 2, eval_col)]
    colnames(best) <- c("query", "target", "evalue")
    return(best)
  }

  fwd <- parse_best(forward_file)
  rev <- parse_best(reverse_file)

  cat(sprintf("Forward best hits: %d queries\n", nrow(fwd)))
  cat(sprintf("Reverse best hits: %d queries\n", nrow(rev)))

  # Find reciprocal pairs: (A→B) and (B→A)
  # Merge on: forward.query == reverse.target AND forward.target == reverse.query
  merged <- merge(fwd, rev, by.x = c("query", "target"), by.y = c("target", "query"))
  rbh_ids <- unique(merged$query)  # query from forward = candidate IDs

  cat(sprintf("Reciprocal best hits: %d pairs\n", length(rbh_ids)))

  # Intersect with candidate list if provided
  if (!is.null(candidates_file) && file.exists(candidates_file)) {
    candidates <- unique(readLines(candidates_file))
    candidates <- candidates[candidates != ""]
    rbh_ids <- intersect(rbh_ids, candidates)
    cat(sprintf("After intersection with candidates: %d\n", length(rbh_ids)))
  }

  # Write output
  writeLines(rbh_ids, output_file)
  cat(sprintf("\nFiltered IDs written to: %s\n", output_file))
  cat(sprintf("Total: %d sequences\n", length(rbh_ids)))
}

# ----- Main -----
if (method == "orthofinder") {
  method_orthofinder()
} else if (method == "rbh") {
  method_rbh()
} else {
  stop(sprintf("Unknown method: %s. Use 'orthofinder' or 'rbh'.", method))
}

cat("\nDone Step 4.\n")
