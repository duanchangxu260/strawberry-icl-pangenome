#!/bin/bash
# =============================================================================
# Step 3: BLASTP against Arabidopsis thaliana proteome
# =============================================================================
# Merges all candidate FASTA files, BLASTs against A. thaliana proteome,
# and extracts best hits per query with optional quality filtering.
#
# Usage:
#   bash 03_blast_arabidopsis.sh <input_dir> <ath_proteome.fa> [options]
#   bash 03_blast_arabidopsis.sh /mnt/e/ogg ath.fa
#   bash 03_blast_arabidopsis.sh /mnt/e/ogg ath.fa --evalue 1e-5 --threads 8
#
# Input:
#   *_ICL_candidates.fasta files in <input_dir>
#   A. thaliana proteome FASTA (e.g., TAIR10_pep_20101214.fasta)
#
# Output:
#   all_candidates.fa          - Merged candidate sequences
#   all_vs_ath.txt             - Raw BLASTP results (outfmt 6)
#   best_hits.txt              - Best hit per query
#   *_ICL_blast_filtered.fasta - Per-species filtered sequences
# =============================================================================

set -euo pipefail

INPUT_DIR="${1:-}"
ATH_FASTA="${2:-}"
EVALUE="${3:-1e-5}"
THREADS="${4:-4}"
MIN_PIDENT="${5:-40}"
MIN_LENGTH="${6:-100}"

if [ -z "$INPUT_DIR" ] || [ -z "$ATH_FASTA" ]; then
    echo "Usage: $0 <input_dir> <ath_proteome.fa> [evalue] [threads] [min_pident] [min_length]"
    echo ""
    echo "Arguments:"
    echo "  input_dir        Directory containing *_ICL_candidates.fasta files"
    echo "  ath_proteome.fa  Arabidopsis thaliana proteome FASTA file"
    echo "  evalue           E-value threshold (default: 1e-5)"
    echo "  threads          Number of CPU threads (default: 4)"
    echo "  min_pident       Minimum percent identity (default: 40)"
    echo "  min_length       Minimum alignment length in aa (default: 100)"
    echo ""
    echo "Example:"
    echo "  $0 /mnt/e/ogg /path/to/TAIR10_pep.fasta 1e-5 8 40 100"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$ATH_FASTA" ]; then
    echo "Error: Arabidopsis proteome FASTA not found: $ATH_FASTA"
    echo ""
    echo "Download TAIR10 proteome:"
    echo "  wget https://www.arabidopsis.org/download_files/Proteins/TAIR10_pep_20101214.fasta -O ath.fa"
    exit 1
fi

# Check for BLAST+
if ! command -v makeblastdb &> /dev/null || ! command -v blastp &> /dev/null; then
    echo "Error: BLAST+ not installed."
    echo "Install with: sudo apt update && sudo apt install ncbi-blast+ -y"
    exit 1
fi

cd "$INPUT_DIR"

ATH_DIR=$(dirname "$ATH_FASTA")
ATH_BASE=$(basename "$ATH_FASTA" .fasta)
ATH_BASE=$(basename "$ATH_BASE" .fa)
ATH_DB="${ATH_DIR}/${ATH_BASE}_db"

echo "=== Step 3: BLASTP against Arabidopsis thaliana ==="
echo "Arabidopsis proteome: $ATH_FASTA"
echo "E-value threshold: $EVALUE"
echo "Min identity: ${MIN_PIDENT}%"
echo "Min length: ${MIN_LENGTH}aa"
echo "Threads: $THREADS"
echo ""

# ----- 3a: Merge all candidate FASTA files -----
echo "--- 3a: Merging candidate FASTA files ---"
find . -name "*_ICL_candidates.fasta" -exec cat {} + > all_candidates.fa
total_seqs=$(grep -c "^>" all_candidates.fa 2>/dev/null || echo 0)
echo "Merged $total_seqs sequences into all_candidates.fa"
echo ""

if [ "$total_seqs" -eq 0 ]; then
    echo "Error: No candidate sequences found."
    exit 1
fi

# ----- 3b: Build Arabidopsis BLAST database -----
echo "--- 3b: Building BLAST database ---"
makeblastdb -in "$ATH_FASTA" -dbtype prot -out "$ATH_DB" 2>&1 | tail -1
echo ""

# ----- 3c: Run BLASTP -----
echo "--- 3c: Running BLASTP ---"
blastp -query all_candidates.fa \
    -db "$ATH_DB" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue "$EVALUE" \
    -num_threads "$THREADS" \
    -out all_vs_ath.txt

hit_count=$(wc -l < all_vs_ath.txt)
echo "BLASTP complete: $hit_count raw hits"
echo ""

if [ "$hit_count" -eq 0 ]; then
    echo "Warning: No BLAST hits found."
    exit 0
fi

# ----- 3d: Extract best hits per query -----
echo "--- 3d: Extracting best hits ---"
# For each query (column 1), keep the line with the lowest E-value (column 11)
awk '!($1 in best) || $11 < eval[$1] { best[$1]=$0; eval[$1]=$11 }
     END { for (q in best) print best[q] }' all_vs_ath.txt | sort > best_hits.txt

best_count=$(wc -l < best_hits.txt)
echo "Best hits: $best_count queries"
echo ""

# ----- 3e: Quality filtering -----
echo "--- 3e: Quality filtering (pident >= ${MIN_PIDENT}, length >= ${MIN_LENGTH}) ---"
awk -v pid="$MIN_PIDENT" -v len="$MIN_LENGTH" \
    '$3 >= pid && $4 >= len {print $0}' best_hits.txt > best_hits_filtered.txt
filtered_count=$(wc -l < best_hits_filtered.txt)
echo "After filtering: $filtered_count hits"
echo ""

# ----- 3f (optional): Per-species BLAST filtering using Arabidopsis ICL -----
# If an Arabidopsis ICL sequence file exists in the PF00463.hmm directory,
# do species-by-species reciprocal BLAST
ATH_ICL_FA=$(find . -path "*/拟南芥/AtICL.fasta" 2>/dev/null | head -1)
if [ -n "$ATH_ICL_FA" ] && [ -f "$ATH_ICL_FA" ]; then
    echo "--- 3f: Per-species BLAST against Arabidopsis ICL ---"
    ATH_ICL_DIR=$(dirname "$ATH_ICL_FA")
    ATH_ICL_DB="${ATH_ICL_DIR}/ath_icl_db"

    # Build Arabidopsis ICL database if not already existing
    if [ ! -f "${ATH_ICL_DB}.pin" ]; then
        makeblastdb -in "$ATH_ICL_FA" -dbtype prot -out "$ATH_ICL_DB" 2>&1 | tail -1
    fi

    find . -name "*_ICL_candidates.fasta" | sort | while read -r query; do
        dir=$(dirname "$query")
        base=$(basename "$query" _ICL_candidates.fasta)
        out="${dir}/${base}_ICL_blast_filtered.fasta"

        echo "  Processing: $(basename "$query")"

        blastp -query "$query" \
            -db "$ATH_ICL_DB" \
            -outfmt "6 qseqid sseqid evalue pident length" \
            -evalue "$EVALUE" \
            -num_threads "$THREADS" > "${dir}/blast_tmp.txt" 2>/dev/null || true

        # Extract passing query IDs
        awk -v ev="$EVALUE" -v pid="$MIN_PIDENT" -v len="$MIN_LENGTH" \
            '$3+0 <= ev+0 && $4+0 >= pid+0 && $5+0 >= len+0 {print $1}' \
            "${dir}/blast_tmp.txt" | sort -u > "${dir}/hit_ids.txt"

        if [ -s "${dir}/hit_ids.txt" ]; then
            seqkit grep -f "${dir}/hit_ids.txt" "$query" -o "$out" 2>/dev/null || \
                awk -v idfile="${dir}/hit_ids.txt" \
                    -f "$(dirname "$(readlink -f "$0")")/extract_fasta_by_ids.awk" \
                    "$query" > "$out"
            kept=$(grep -c "^>" "$out" 2>/dev/null || echo 0)
            echo "    → $kept sequences kept"
        else
            echo "    → No hits passed filter"
            touch "$out"
        fi

        rm -f "${dir}/blast_tmp.txt" "${dir}/hit_ids.txt"
    done
    echo ""
fi

echo "Done Step 3."
echo ""
echo "Output files:"
echo "  all_candidates.fa           - Merged candidate sequences ($total_seqs seqs)"
echo "  all_vs_ath.txt              - Raw BLASTP results ($hit_count hits)"
echo "  best_hits.txt               - Best BLAST hit per query ($best_count queries)"
echo "  best_hits_filtered.txt      - Quality-filtered best hits ($filtered_count queries)"
