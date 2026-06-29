#!/bin/bash
# =============================================================================
# Step 1: Filter HMMER hmmsearch results by E-value (Shell/AWK version)
# =============================================================================
# Parses HMMER --tblout output (*.ICL.txt) and extracts sequence IDs
# with E-value <= threshold (default 1e-5).
#
# Usage:
#   bash 01_filter_hmm.sh <input_dir> [evalue_threshold]
#   bash 01_filter_hmm.sh /mnt/e/ogg 1e-5
#   bash 01_filter_hmm.sh /mnt/e/ogg 1e-10
#
# Input:
#   *.ICL.txt files in <input_dir> (HMMER --tblout format)
#
# Output:
#   *_filtered_ids.txt files in the same directories as inputs
# =============================================================================

set -euo pipefail

# Parse arguments
INPUT_DIR="${1:-}"
EVALUE_THRESHOLD="${2:-1e-5}"

if [ -z "$INPUT_DIR" ]; then
    echo "Usage: $0 <input_dir> [evalue_threshold]"
    echo "Example: $0 /mnt/e/ogg 1e-5"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory not found: $INPUT_DIR"
    exit 1
fi

cd "$INPUT_DIR"

echo "=== Step 1: Filter HMMER results by E-value ==="
echo "Input directory: $INPUT_DIR"
echo "E-value threshold: $EVALUE_THRESHOLD"
echo ""

# Count files
file_count=$(find . -name "*.ICL.txt" | wc -l)
echo "Found $file_count *.ICL.txt file(s)"
echo ""

if [ "$file_count" -eq 0 ]; then
    echo "Warning: No *.ICL.txt files found. Exiting."
    exit 0
fi

# Process each file
# HMMER tblout format (0-indexed columns):
#   [0] target name
#   [4] full sequence E-value
#   [5] full sequence score
#
# Skipping comment lines (starting with #), filtering where col[4] <= threshold

total_filtered=0

find . -name "*.ICL.txt" | sort | while read -r f; do
    dir=$(dirname "$f")
    base=$(basename "$f" .ICL.txt)
    out="${dir}/${base}_filtered_ids.txt"

    # AWK: skip comments (#), check column 5 (1-indexed = $5) for E-value
    # $5 corresponds to column index 4 (full sequence E-value)
    awk -v ev="$EVALUE_THRESHOLD" '
        !/^#/ {
            if (NF >= 5 && $5+0 <= ev+0) {
                print $1
            }
        }
    ' "$f" | sort -u > "$out"

    count=$(wc -l < "$out")
    echo "  $(basename "$f"): $count IDs passed filter → $(basename "$out")"
    total_filtered=$((total_filtered + count))
done

echo ""
echo "Done Step 1."
