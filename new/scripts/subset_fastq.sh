#!/bin/bash

# Bash strict mode
set -euo pipefail

#
# Description:
# This script filters gzipped R1 and R2 FASTQ files based on a list of read names.
# It takes a two-column mapping file (readname, genename) and two FASTQ files as input.
# It outputs a single-column file of gene names and two new, filtered FASTQ files.
#
# Dependencies:
#   - seqtk
#   - gzip
#
# Usage:
#   ./filter_fastq_by_name.sh <read_gene_map.txt> <reads_R1.fastq.gz> <reads_R2.fastq.gz>
#
# Example:
#   ./filter_fastq_by_name.sh my_reads_to_genes.txt input_1.fastq.gz input_2.fastq.gz
#

# --- Functions ---

# Function to display usage information
usage() {
    echo "Usage: $0 <read_gene_map_file> <r1_fastq_gz> <r2_fastq_gz>"
    echo ""
    echo "Arguments:"
    echo "  read_gene_map_file    A tab or space-separated file with two columns: read names and gene names."
    echo "  r1_fastq_gz           The gzipped R1 FASTQ file."
    echo "  r2_fastq_gz           The gzipped R2 FASTQ file."
    exit 1
}

# Function to check for required command-line tools
check_dependencies() {
    for cmd in seqtk gzip; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: Required command '$cmd' is not found."
            echo "Please install it and ensure it's in your system's PATH."
            exit 1
        fi
    done
}


# --- Main Script ---

# 1. Argument Validation
if [ "$#" -ne 3 ]; then
    echo "Error: Incorrect number of arguments."
    usage
fi

# 2. Check for dependencies
check_dependencies

# 3. Assign arguments to variables for clarity
READ_GENE_MAP="$1"
R1_FASTQ="$2"
R2_FASTQ="$3"

# 4. Check if input files exist
for f in "$READ_GENE_MAP" "$R1_FASTQ" "$R2_FASTQ"; do
    if [ ! -f "$f" ]; then
        echo "Error: Input file not found: $f"
        exit 1
    fi
done

# 5. Define output file names
GENE_NAMES_OUTPUT="filtered_gene_names.txt"
R1_FASTQ_OUTPUT="filtered_R1.fastq.gz"
R2_FASTQ_OUTPUT="filtered_R2.fastq.gz"
READ_NAMES_TMP=$(mktemp) # Create a secure temporary file for read names

echo "--> Script started."
echo "--> Read/Gene Map: $READ_GENE_MAP"
echo "--> R1 FASTQ:      $R1_FASTQ"
echo "--> R2 FASTQ:      $R2_FASTQ"
echo ""

# 6. Extract read names (column 1) and gene names (column 2)
echo "Step 1: Extracting read and gene names..."
# awk is used for flexibility with tab or space delimiters
awk '{print $1}' "$READ_GENE_MAP" > "$READ_NAMES_TMP"
awk '{print $2}' "$READ_GENE_MAP" > "$GENE_NAMES_OUTPUT"
echo "Done. Gene names saved to '$GENE_NAMES_OUTPUT'"

# 7. Filter the R1 and R2 FASTQ files using seqtk
echo "Step 2: Filtering R1 FASTQ file with seqtk..."
seqtk subseq "$R1_FASTQ" "$READ_NAMES_TMP" | gzip > "$R1_FASTQ_OUTPUT"
echo "Done. Filtered R1 saved to '$R1_FASTQ_OUTPUT'"

echo "Step 3: Filtering R2 FASTQ file with seqtk..."
seqtk subseq "$R2_FASTQ" "$READ_NAMES_TMP" | gzip > "$R2_FASTQ_OUTPUT"
echo "Done. Filtered R2 saved to '$R2_FASTQ_OUTPUT'"

# 8. Clean up the temporary file
rm "$READ_NAMES_TMP"
echo ""
echo "--> All tasks completed successfully!"
echo "--> Output files created:"
echo "    1. $GENE_NAMES_OUTPUT"
echo "    2. $R1_FASTQ_OUTPUT"
echo "    3. $R2_FASTQ_OUTPUT"
