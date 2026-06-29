#!/bin/bash
# =============================================================================
# ICL Pan-Genome Analysis — Complete Pipeline
# =============================================================================
# One-click execution of all four steps for ICL gene family identification.
#
# Usage:
#   bash run_pipeline.sh <input_dir> <ath_proteome.fa> [options]
#
#   bash run_pipeline.sh /mnt/e/ogg ath.fa
#   bash run_pipeline.sh /mnt/e/ogg ath.fa --evalue 1e-5 --threads 8
#   bash run_pipeline.sh /mnt/e/ogg ath.fa --skip-ogg  # Skip OGG filtering
#
# Required:
#   <input_dir>        Directory with *.ICL.txt and proteome FASTA files
#   <ath_proteome.fa>  Arabidopsis thaliana proteome FASTA
#
# Options:
#   --evalue FLOAT     E-value threshold for HMMER and BLAST (default: 1e-5)
#   --threads INT      Number of CPU threads (default: 4)
#   --min-pident INT   Minimum BLAST identity % (default: 40)
#   --min-length INT   Minimum BLAST alignment length in aa (default: 100)
#   --skip-blast       Skip BLAST step (Step 3)
#   --skip-ogg         Skip OGG filtering step (Step 4)
#   --orthogroups FILE Path to OrthoFinder Orthogroups.tsv for OGG filtering
#   --dry-run          Print steps without executing
#   -h, --help         Show this help message
# =============================================================================

set -euo pipefail

# ----- Configuration -----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default parameters
EVALUE="1e-5"
THREADS="4"
MIN_PIDENT="40"
MIN_LENGTH="100"
SKIP_BLAST=false
SKIP_OGG=false
ORTHOGROUPS=""
DRY_RUN=false
ATH_FASTA=""
INPUT_DIR=""

# ----- Parse arguments -----
while [[ $# -gt 0 ]]; do
    case $1 in
        --evalue)
            EVALUE="$2"; shift 2 ;;
        --threads)
            THREADS="$2"; shift 2 ;;
        --min-pident)
            MIN_PIDENT="$2"; shift 2 ;;
        --min-length)
            MIN_LENGTH="$2"; shift 2 ;;
        --skip-blast)
            SKIP_BLAST=true; shift ;;
        --skip-ogg)
            SKIP_OGG=true; shift ;;
        --orthogroups)
            ORTHOGROUPS="$2"; shift 2 ;;
        --dry-run)
            DRY_RUN=true; shift ;;
        -h|--help)
            head -30 "$0" | grep -A 100 "^#" | sed 's/^# \?//'
            exit 0 ;;
        -*)
            echo "Unknown option: $1" >&2; exit 1 ;;
        *)
            if [ -z "$INPUT_DIR" ]; then
                INPUT_DIR="$1"
            elif [ -z "$ATH_FASTA" ]; then
                ATH_FASTA="$1"
            else
                echo "Too many positional arguments" >&2
                exit 1
            fi
            shift ;;
    esac
done

# ----- Validation -----
if [ -z "$INPUT_DIR" ] || [ -z "$ATH_FASTA" ]; then
    echo "Usage: $0 <input_dir> <ath_proteome.fa> [options]"
    echo "Run '$0 --help' for details."
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$ATH_FASTA" ] && [ "$SKIP_BLAST" = false ]; then
    echo "Error: Arabidopsis proteome FASTA not found: $ATH_FASTA"
    exit 1
fi

# ----- Print configuration -----
echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║     ICL Pan-Genome Analysis Pipeline                         ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Input directory:     $INPUT_DIR"
echo "A. thaliana FASTA:   $ATH_FASTA"
echo "E-value threshold:   $EVALUE"
echo "BLAST threads:       $THREADS"
echo "Min identity:        ${MIN_PIDENT}%"
echo "Min align length:    ${MIN_LENGTH}aa"
echo "Skip BLAST:          $SKIP_BLAST"
echo "Skip OGG filter:     $SKIP_OGG"
echo ""

if [ "$DRY_RUN" = true ]; then
    echo "[DRY RUN MODE] No commands will be executed."
    echo ""
    echo "Would run:"
    echo "  1. bash ${SCRIPT_DIR}/01_filter_hmm.sh '$INPUT_DIR' '$EVALUE'"
    echo "  2. bash ${SCRIPT_DIR}/02_extract_sequences.sh '$INPUT_DIR'"
    if [ "$SKIP_BLAST" = false ]; then
        echo "  3. bash ${SCRIPT_DIR}/03_blast_arabidopsis.sh '$INPUT_DIR' '$ATH_FASTA' '$EVALUE' '$THREADS' '$MIN_PIDENT' '$MIN_LENGTH'"
    fi
    if [ "$SKIP_OGG" = false ] && [ -n "$ORTHOGROUPS" ]; then
        echo "  4. Rscript ${SCRIPT_DIR}/04_ogg_filter.R orthofinder --orthogroups '$ORTHOGROUPS' --reference AT3G21720 --candidates best_hits.txt"
    fi
    exit 0
fi

# ----- Step 1: E-value filtering -----
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  Step 1/4: E-value filtering                                 ║"
echo "╚══════════════════════════════════════════════════════════════╝"
bash "${SCRIPT_DIR}/01_filter_hmm.sh" "$INPUT_DIR" "$EVALUE"

# ----- Step 2: Sequence extraction -----
echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  Step 2/4: Protein sequence extraction                       ║"
echo "╚══════════════════════════════════════════════════════════════╝"
bash "${SCRIPT_DIR}/02_extract_sequences.sh" "$INPUT_DIR"

# ----- Step 3: BLAST against Arabidopsis -----
if [ "$SKIP_BLAST" = false ]; then
    echo ""
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║  Step 3/4: BLASTP against Arabidopsis thaliana               ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    bash "${SCRIPT_DIR}/03_blast_arabidopsis.sh" "$INPUT_DIR" "$ATH_FASTA" "$EVALUE" "$THREADS" "$MIN_PIDENT" "$MIN_LENGTH"
else
    echo ""
    echo "Step 3/4: BLAST — SKIPPED"
fi

# ----- Step 4: OGG Filtering -----
if [ "$SKIP_OGG" = false ]; then
    echo ""
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║  Step 4/4: Orthogroup (OGG) filtering                        ║"
    echo "╚══════════════════════════════════════════════════════════════╝"

    if [ -n "$ORTHOGROUPS" ] && [ -f "$ORTHOGROUPS" ]; then
        # Method A: OrthoFinder
        echo "Using OrthoFinder orthogroups: $ORTHOGROUPS"
        Rscript "${SCRIPT_DIR}/04_ogg_filter.R" orthofinder \
            --orthogroups "$ORTHOGROUPS" \
            --reference "AT3G21720" \
            --candidates "${INPUT_DIR}/best_hits.txt" \
            --output "${INPUT_DIR}/ogg_filtered_ids.txt"
    elif [ -f "${INPUT_DIR}/all_vs_ath.txt" ]; then
        # Method B: RBH (if reverse BLAST was done)
        echo "Note: No orthogroups file provided. To complete OGG filtering:"
        echo "  1. Run eggNOG-mapper on all_candidates.fa"
        echo "  2. Or provide --orthogroups path to Orthogroups.tsv"
        echo "  3. Or use RBH method with reverse BLAST results"
        echo ""
        echo "Skipping OGG filtering for now."
    else
        echo "Warning: Cannot perform OGG filtering without orthogroups or BLAST results."
        echo "Please provide --orthogroups, or run eggNOG-mapper on all_candidates.fa"
    fi
else
    echo ""
    echo "Step 4/4: OGG filtering — SKIPPED"
fi

# ----- Summary -----
echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  Pipeline Complete                                           ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Output files in $INPUT_DIR:"
echo ""
echo "  Step 1: *_filtered_ids.txt"
find "$INPUT_DIR" -name "*_filtered_ids.txt" | sort | while read -r f; do
    count=$(wc -l < "$f")
    echo "    $(basename "$f") ($count IDs)"
done

echo ""
echo "  Step 2: *_ICL_candidates.fasta"
find "$INPUT_DIR" -name "*_ICL_candidates.fasta" | sort | while read -r f; do
    count=$(grep -c "^>" "$f" 2>/dev/null || echo 0)
    echo "    $(basename "$f") ($count sequences)"
done

if [ -f "${INPUT_DIR}/all_candidates.fa" ]; then
    echo ""
    echo "  Step 3:"
    echo "    all_candidates.fa ($(grep -c '^>' ${INPUT_DIR}/all_candidates.fa) sequences)"
    [ -f "${INPUT_DIR}/best_hits.txt" ] && echo "    best_hits.txt ($(wc -l < ${INPUT_DIR}/best_hits.txt) hits)"
fi

if [ -f "${INPUT_DIR}/ogg_filtered_ids.txt" ]; then
    echo ""
    echo "  Step 4:"
    echo "    ogg_filtered_ids.txt ($(wc -l < ${INPUT_DIR}/ogg_filtered_ids.txt) IDs)"
fi

echo ""
echo "Done!"
