#!/usr/bin/awk -f
# =============================================================================
# AWK helper: Extract FASTA sequences by ID list
# =============================================================================
# Reads a list of target IDs from a file and extracts matching sequences
# from a FASTA file on stdin.
#
# Usage:
#   awk -v idfile="filtered_ids.txt" -f extract_fasta_by_ids.awk input.fa > output.fasta
#
# The ID file should contain one sequence ID per line.
# Matching is done against the first space-delimited token in the FASTA header.
# =============================================================================

BEGIN {
    # Load target IDs into a hash
    while ((getline < idfile) > 0) {
        ids[$1] = 1
    }
    # Close the idfile to avoid conflicts
    close(idfile)
}

/^>/ {
    p = 0  # Reset print flag
    # Extract the first token from the header (remove the '>')
    header = substr($0, 2)
    split(header, parts, " ")
    seq_id = parts[1]

    if (seq_id in ids) {
        p = 1
    }
}

# Print lines when p flag is set
p { print }
