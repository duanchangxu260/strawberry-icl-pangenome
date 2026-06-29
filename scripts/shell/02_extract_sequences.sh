#!/bin/bash
# =============================================================================
# Step 2: Extract protein sequences from FASTA using filtered ID lists
# =============================================================================
# Reads *_filtered_ids.txt files and extracts corresponding sequences from
# species proteome FASTA files. Uses seqkit if available, falls back to AWK.
#
# Usage:
#   bash 02_extract_sequences.sh <input_dir>
#   bash 02_extract_sequences.sh /mnt/e/ogg
#
# Input:
#   *_filtered_ids.txt (from Step 1)
#   Species proteome FASTA files (*.pep.fasta, *.fa, *.fasta, *.faa, *.pep)
#
# Output:
#   *_ICL_candidates.fasta in the same directories
# =============================================================================

set -euo pipefail

INPUT_DIR="${1:-}"

if [ -z "$INPUT_DIR" ]; then
    echo "Usage: $0 <input_dir>"
    echo "Example: $0 /mnt/e/ogg"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory not found: $INPUT_DIR"
    exit 1
fi

cd "$INPUT_DIR"

# Check if seqkit is available
if command -v seqkit &> /dev/null; then
    EXTRACTOR="seqkit"
    echo "Using seqkit for sequence extraction"
else
    EXTRACTOR="awk"
    echo "seqkit not found, using AWK (slower but no installation needed)"
    echo "Tip: Install seqkit for faster extraction:"
    echo "  conda install -c bioconda seqkit"
    echo ""
fi

echo "=== Step 2: Extract protein sequences ==="
echo ""

# Find all ID list files
id_count=$(find . -name "*_filtered_ids.txt" | wc -l)
echo "Found $id_count *_filtered_ids.txt file(s)"
echo ""

total_extracted=0

find . -name "*_filtered_ids.txt" | sort | while read -r idfile; do
    dir=$(dirname "$idfile")
    base=$(basename "$idfile" _filtered_ids.txt)
    out="${dir}/${base}_ICL_candidates.fasta"

    # Find FASTA file in the same directory
    # Priority: .pep.fasta > .fa > .fasta > .faa
    fa=""
    for ext in ".pep.fasta" ".fa" ".fasta" ".faa" ".pep" ".protein.fasta"; do
        fa=$(ls "$dir"/*"$ext" 2>/dev/null | head -1) || true
        [ -n "$fa" ] && break
    done

    if [ -z "$fa" ]; then
        echo "  WARNING: No FASTA file found in $dir, skipping $(basename "$idfile")"
        continue
    fi

    # Extract sequences
    if [ "$EXTRACTOR" = "seqkit" ]; then
        # -r: treat IDs as regex, -i: case-insensitive, -f: read IDs from file
        seqkit grep -r -i -f "$idfile" "$fa" -o "$out" 2>/dev/null || {
            # If regex mode fails, try exact match (without -r)
            seqkit grep -i -f "$idfile" "$fa" -o "$out" 2>/dev/null
        }
    else
        # Use AWK helper script
        script_dir=$(dirname "$(readlink -f "$0")")
        awk -v idfile="$idfile" -f "${script_dir}/extract_fasta_by_ids.awk" "$fa" > "$out"
    fi

    count=$(grep -c "^>" "$out" 2>/dev/null || echo 0)
    echo "  $(basename "$idfile"): $count sequences → $(basename "$out")"
    total_extracted=$((total_extracted + count))
done

echo ""
echo "Done Step 2."
