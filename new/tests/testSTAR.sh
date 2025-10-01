#!/bin/bash

# Run STAR with CB/UB tag table export
# Based on production parameters from SC2300771

set -euo pipefail

# Configuration
NEW_STAR_BINARY="./bin/Linux_x86_64/STAR"                     # Newly built binary in current directory
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
THREADS=24
OUTPUT_BASE="./testing"

NEW_DIR="${OUTPUT_BASE}"
TEMP_DIR="/storage/tmp/testing"

# Parse command line arguments for extra flags
EXTRA_FLAGS=()
while [[ $# -gt 0 ]]; do
    EXTRA_FLAGS+=("$1")
    shift
done

#clean up the temp directory
rm -rf $TEMP_DIR


# Set ulimit for STAR
ulimit -n 100000 || true

# Helper: clean dir
cleanup_dir() {
    local dir="$1"
    echo "Cleaning up: $dir"
    rm -rf "$dir" "${dir}_STARtmp" 2>/dev/null || true
    mkdir -p "$dir"
}

# Check requirements
echo "=== Checking inputs/binaries ==="
[[ -f "./filtered_R1.fastq.gz" ]] || { echo "ERROR: R1 file not found: ./filtered_R1.fastq.gz"; exit 1; }
[[ -f "./filtered_R2.fastq.gz" ]] || { echo "ERROR: R2 file not found: ./filtered_R2.fastq.gz"; exit 1; }
[[ -f "$WHITELIST" ]] || { echo "ERROR: Whitelist not found: $WHITELIST"; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: Genome dir not found: $GENOME_DIR"; exit 1; }
[[ -f "$NEW_STAR_BINARY" ]] || { echo "ERROR: New STAR binary not found: $NEW_STAR_BINARY"; exit 1; }

# Set R2/R1 file paths
R2_FILES="./filtered_R2.fastq.gz"
R1_FILES="./filtered_R1.fastq.gz"

echo "R2: $R2_FILES"
echo "R1: $R1_FILES"

# Common STAR params
COMMON_PARAMS=(
    --runThreadN $THREADS
    --outTmpDir $TEMP_DIR
    --quantMode GeneCounts
    --soloType CB_UMI_Simple
    --soloCBlen 16
    --soloUMIlen 12
    --soloUMIstart 17
    --soloCBstart 1
    --soloCBwhitelist "$WHITELIST"
    --genomeDir "$GENOME_DIR"
    --limitIObufferSize 50000000 50000000
    --outSJtype None
    --outBAMcompression 6
    --soloMultiMappers Rescue
    --alignIntronMax 1
    --alignMatesGapMax 0
    --outFilterMismatchNmax 6
    --outFilterMismatchNoverReadLmax 1.0
    --outFilterMatchNmin 25
    --soloBarcodeReadLength 0
    --outSAMunmapped None
    --outFilterMatchNminOverLread 0
    --outFilterMultimapNmax 10000
    --outFilterMultimapScoreRange 1
    --outSAMmultNmax 10000
    --winAnchorMultimapNmax 200
    --outSAMprimaryFlag AllBestScore
    --outFilterScoreMin 0
    --outFilterScoreMinOverLread 0
    --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN gx gn ZG ZX
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloCellFilter None
    --clipAdapterType CellRanger4
    --soloFeatures Gene GeneFull
    --soloStrand Unstranded
    --readFilesIn "$R2_FILES" "$R1_FILES"
    --alignEndsType Local
    --readFilesCommand zcat
)

NEW_TAG_TABLE_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
    --soloAddTagsToUnsorted yes
    --soloWriteTagTable Default
    "${EXTRA_FLAGS[@]}"
)

# Run new STAR with tag table export
cleanup_dir "$NEW_DIR"

echo "=== Running new STAR with CB/UB tag table export ==="
if [[ ${#EXTRA_FLAGS[@]} -gt 0 ]]; then
    echo "Extra flags: ${EXTRA_FLAGS[*]}"
fi
echo "$NEW_STAR_BINARY \
    "${NEW_TAG_TABLE_PARAMS[@]}" \
    --outFileNamePrefix ${NEW_DIR}/"
eval  "$NEW_STAR_BINARY \
    "${NEW_TAG_TABLE_PARAMS[@]}" \
    --outFileNamePrefix ${NEW_DIR}/"

echo "Run finished"
