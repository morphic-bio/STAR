#!/bin/bash

# Run STAR with CB/UB tag table export
# Based on production parameters from SC2300771
#
# =============================================================================
# ZG/ZX BAM Tags Documentation
# =============================================================================
#
# This script includes support for custom ZG (gene set) and ZX (overlap status) 
# BAM tags that provide detailed gene annotation information for each read.
#
# ZG Tag (Gene Set):
# ------------------
# - Format: ZG:Z:<gene_list> where <gene_list> is comma-separated Ensembl gene IDs
# - Example: ZG:Z:ENSG00000103024,ENSG00000171824
# - Empty value: ZG:Z:- (when no genes overlap the read)
# - Contains all genes that overlap with the read based on genomic coordinates
# - Uses GeneFull annotation mode for comprehensive gene detection
#
# ZX Tag (Overlap Status):
# ------------------------
# - Format: ZX:Z:<overlap_type> where <overlap_type> describes the genomic context
# - Valid values:
#   * "none"      - Read doesn't overlap any annotated features (intergenic)
#   * "exonic"    - Read overlaps with exonic regions (includes antisense)
#   * "intronic"  - Read overlaps with intronic regions (includes antisense)  
#   * "spanning"  - Read spans multiple feature types or unknown overlap
# - Provides genomic context information for downstream analysis
#
# Implementation Details:
# -----------------------
# - ZG/ZX tags are BAM-only (not emitted in SAM output due to complexity)
# - Requires --soloFeatures to include "GeneFull" for proper gene detection
# - Recommends --soloStrand Unstranded to avoid strand specificity issues
# - Tags are populated using existing STAR gene annotation infrastructure
# - Minimal performance impact as they reuse existing annotation data structures
#
# Usage Requirements:
# -------------------
# 1. Add "ZG ZX" to --outSAMattributes parameter (already included below)
# 2. Ensure --soloFeatures includes "GeneFull" (already configured)
# 3. Use BAM output format (--outSAMtype BAM)
# 4. Gene annotation must be available in genome index
#
# Example Output:
# ---------------
# Read mapping to NME3 gene:     ZG:Z:ENSG00000103024  ZX:Z:exonic
# Read mapping to EXOSC10 gene:  ZG:Z:ENSG00000171824  ZX:Z:intronic  
# Read with no gene overlap:     ZG:Z:-                ZX:Z:none
# Read mapping multiple genes:   ZG:Z:ENSG00001,ENSG00002  ZX:Z:exonic
#
# Troubleshooting:
# ----------------
# - Empty ZG tags (ZG:Z:-): Check gene annotation coverage and strand settings
# - Missing ZG/ZX tags: Ensure "ZG ZX" is in --outSAMattributes
# - Compilation errors: Run 'make clean && make' to rebuild with ZG/ZX support
#
# For validation, use: python3 validate_zg_zx.py output.bam allowed_genes.txt
# =============================================================================

set -euo pipefail

# Configuration
NEW_STAR_BINARY="./bin/Linux_x86_64/STAR"                     # Newly built binary in current directory
SAMPLE_ID="SC2300771"
BASE_DIR="/storage/downsampled/${SAMPLE_ID}"
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
THREADS=24
OUTPUT_BASE="/storage/SC2300771_subset"

NEW_DIR="${OUTPUT_BASE}/${SAMPLE_ID}"
TEMP_DIR="/storage/tmp/${SAMPLE_ID}"

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
[[ -d "$BASE_DIR" ]] || { echo "ERROR: Input dir not found: $BASE_DIR"; exit 1; }
[[ -f "$WHITELIST" ]] || { echo "ERROR: Whitelist not found: $WHITELIST"; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: Genome dir not found: $GENOME_DIR"; exit 1; }
[[ -f "$NEW_STAR_BINARY" ]] || { echo "ERROR: New STAR binary not found: $NEW_STAR_BINARY"; exit 1; }

# Build R2/R1 file lists (8 lanes)
R2_FILES=""; R1_FILES=""
for lane in L001 L002 L003 L004 L005 L006 L007 L008; do
    r2_glob="${BASE_DIR}/${SAMPLE_ID}_*_${lane}_R2_001.fastq.gz"
    r1_glob="${BASE_DIR}/${SAMPLE_ID}_*_${lane}_R1_001.fastq.gz"
    if ls $r2_glob 1>/dev/null 2>&1; then
        [[ -n "$R2_FILES" ]] && R2_FILES+="",
        R2_FILES+="$(ls $r2_glob)"
    fi
    if ls $r1_glob 1>/dev/null 2>&1; then
        [[ -n "$R1_FILES" ]] && R1_FILES+="",
        R1_FILES+="$(ls $r1_glob)"
    fi
done
[[ -n "$R2_FILES" && -n "$R1_FILES" ]] || { echo "ERROR: Could not find R1/R2 FASTQs in $BASE_DIR"; exit 1; }

echo "R2: $R2_FILES"
echo "R1: $R1_FILES"

# Common STAR params
# Note: ZG/ZX BAM tags are configured via the following parameters:
# - --outSAMattributes includes "ZG ZX" for gene set and overlap status tags
# - --soloFeatures includes "GeneFull" for comprehensive gene annotation
# - --soloStrand Unstranded recommended to avoid strand specificity issues
# - --outSAMtype BAM required (ZG/ZX are BAM-only tags)
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
    --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN gx gn ZG ZX  # ZG=gene_set, ZX=overlap_status
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloCellFilter None
    --clipAdapterType CellRanger4
    --soloFeatures Gene GeneFull                                     # GeneFull required for ZG/ZX tags
    --soloStrand Unstranded                                          # Recommended for ZG/ZX tags
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

# =============================================================================
# ZG/ZX Output Interpretation Guide
# =============================================================================
#
# After STAR completes, the output BAM file will contain ZG and ZX tags for each read.
# Use samtools or other BAM processing tools to examine the tags:
#
# Example commands to inspect ZG/ZX tags:
# samtools view output.bam | grep -E "ZG:Z:|ZX:Z:" | head -10
# samtools view output.bam | awk '{for(i=1;i<=NF;i++) if($i~/^ZG:Z:/ || $i~/^ZX:Z:/) print $1"\t"$i}' | head -10
#
# Interpretation examples:
# - ZG:Z:ENSG00000103024 ZX:Z:exonic    → Read maps to NME3 gene, overlaps exonic region
# - ZG:Z:- ZX:Z:none                    → Read doesn't overlap any annotated genes
# - ZG:Z:ENSG00001,ENSG00002 ZX:Z:exonic → Read overlaps multiple genes in exonic regions
#
# For downstream analysis:
# - Use ZG tags to identify gene-specific reads for expression quantification
# - Use ZX tags to filter reads by genomic context (exonic vs intronic vs intergenic)
# - Combine with existing tags (GX, GN) for comprehensive gene annotation analysis
# =============================================================================
