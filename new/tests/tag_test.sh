#!/bin/bash

# Simplified Tag Test - Final validation of STARsolo two-pass unsorted CB/UB tag injection
# This version does NOT re-run the baseline. Baseline must already exist.
# It runs the new STAR with --soloAddTagsToUnsorted and compares to the baseline.
# Based on production parameters from SC2300771

set -euo pipefail

# Configuration
ORIGINAL_STAR_BINARY="/usr/local/bin/STAR"   # Not used here (baseline is not re-run)
NEW_STAR_BINARY="./STAR"                     # Newly built binary in current directory
SAMPLE_ID="SC2300771"
BASE_DIR="/storage/downsampled/${SAMPLE_ID}"
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
THREADS=24
OUTPUT_BASE="$(pwd)/tag_test_output"

# Optional flags
BASELINE_DIR="${OUTPUT_BASE}/baseline"       # Must contain baseline from original STAR
NEW_DIR="${OUTPUT_BASE}/new_with_tags"

# Set ulimit for STAR
ulimit -n 100000 || true

# Helper: clean dir
cleanup_dir() {
    local dir="$1"
    echo "Cleaning up: $dir"
    rm -rf "$dir" "${dir}_STARtmp" 2>/dev/null || true
    mkdir -p "$dir"
}

# Helper: compare BAMs by record count (order may differ)
compare_bams_counts() {
    local bam1="$1"
    local bam2="$2"
    local name1="$3"
    local name2="$4"
    
    echo "Comparing BAM record counts:"
    echo "  $name1: $(du -h "$bam1" | cut -f1)"
    echo "  $name2: $(du -h "$bam2" | cut -f1)"
    
    local count1=$(samtools view -c "$bam1" 2>/dev/null || echo "Error")
    local count2=$(samtools view -c "$bam2" 2>/dev/null || echo "Error")
    echo "  $name1 records: $count1"
    echo "  $name2 records: $count2"
    
    if [[ "$count1" == "$count2" && "$count1" != "Error" ]]; then
        echo "  ✓ Record counts match"
        return 0
    else
        echo "  ✗ Record counts differ"
        return 1
    fi
}

# Helper: compare Solo matrices byte-identically
compare_solo() {
    local solo1="$1"
    local solo2="$2"
    
    local features1="${solo1}/Gene/raw/features.tsv"
    local features2="${solo2}/Gene/raw/features.tsv"
    local barcodes1="${solo1}/Gene/raw/barcodes.tsv"
    local barcodes2="${solo2}/Gene/raw/barcodes.tsv"
    local matrix1="${solo1}/Gene/raw/matrix.mtx"
    local matrix2="${solo2}/Gene/raw/matrix.mtx"

    local ok=true
    if [[ -f "$features1" && -f "$features2" ]] && diff -q "$features1" "$features2" >/dev/null; then
        echo "  ✓ Features.tsv identical"
    else
        echo "  ✗ Features.tsv differ or missing"; ok=false
    fi
    if [[ -f "$barcodes1" && -f "$barcodes2" ]] && diff -q "$barcodes1" "$barcodes2" >/dev/null; then
        echo "  ✓ Barcodes.tsv identical"
    else
        echo "  ✗ Barcodes.tsv differ or missing"; ok=false
    fi
    if [[ -f "$matrix1" && -f "$matrix2" ]] && diff -q "$matrix1" "$matrix2" >/dev/null; then
        echo "  ✓ Matrix.mtx identical"
    else
        echo "  ✗ Matrix.mtx differ or missing"; ok=false
    fi

    $ok && return 0 || return 1
}

# Check requirements
echo "=== Checking inputs/binaries ==="
[[ -d "$BASE_DIR" ]] || { echo "ERROR: Input dir not found: $BASE_DIR"; exit 1; }
[[ -f "$WHITELIST" ]] || { echo "ERROR: Whitelist not found: $WHITELIST"; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: Genome dir not found: $GENOME_DIR"; exit 1; }
[[ -f "$NEW_STAR_BINARY" ]] || { echo "ERROR: New STAR binary not found: $NEW_STAR_BINARY"; exit 1; }

# Baseline must exist
BASELINE_BAM="${BASELINE_DIR}/Aligned.sortedByCoord.out.bam"
BASELINE_SOLO="${BASELINE_DIR}/Solo.out"
if [[ ! -f "$BASELINE_BAM" || ! -d "$BASELINE_SOLO" ]]; then
    echo "ERROR: Baseline missing. Expected:"
    echo "  $BASELINE_BAM"
    echo "  $BASELINE_SOLO"
    echo "Please run the baseline once (e.g., via integration_test.sh Test 1a) and re-run this script."
    exit 1
fi

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
COMMON_PARAMS=(
    --runThreadN $THREADS
    --soloType CB_UMI_Simple
    --soloCBlen 16
    --soloUMIlen 12
    --soloUMIstart 17
    --soloCBstart 1
    --soloBarcodeReadLength 0
    --soloCBwhitelist "$WHITELIST"
    --genomeDir "$GENOME_DIR"
    --limitIObufferSize 50000000 50000000
    --outSJtype None
    --outBAMcompression 6
    --soloMultiMappers Unique
    --alignIntronMax 1
    --alignMatesGapMax 0
    --outFilterMismatchNmax 10
    --outFilterMismatchNoverReadLmax 1.0
    --outFilterMatchNmin 16
    --outSAMunmapped None
    --outFilterMatchNminOverLread 0
    --outFilterMultimapNmax 10000
    --outFilterMultimapScoreRange 1000
    --outSAMmultNmax 10000
    --winAnchorMultimapNmax 200
    --outSAMprimaryFlag AllBestScore
    --outFilterScoreMin 0
    --outFilterScoreMinOverLread 0
    --outSAMattributes NH HI AS nM NM CB UB CR CY UR UY GX GN
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloCellFilter None
    --clipAdapterType CellRanger4
    --soloFeatures Gene
    --readFilesIn "$R2_FILES" "$R1_FILES"
    --alignEndsType Local
    --readFilesCommand zcat
)

NEW_UNSORTED_TAGS_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
    --soloAddTagsToUnsorted yes
)

# Run new STAR with tags
cleanup_dir "$NEW_DIR"

echo "=== Running new STAR with CB/UB tags ==="
$NEW_STAR_BINARY \
    "${NEW_UNSORTED_TAGS_PARAMS[@]}" \
    --outFileNamePrefix "${NEW_DIR}/"

echo "Run finished"

# Validate outputs
NEW_BAM="${NEW_DIR}/Aligned.out.bam"
NEW_SOLO="${NEW_DIR}/Solo.out"

if [[ ! -f "$NEW_BAM" ]]; then
    echo "❌ ERROR: New BAM not created: $NEW_BAM"
    echo "Searching under ${NEW_DIR} for misplaced BAM:"
    find "$NEW_DIR" -maxdepth 2 -type f -name 'Aligned*.bam' -printf '  %p\n' || true
    exit 1
fi
[[ -d "$NEW_SOLO" ]] || { echo "❌ ERROR: New Solo output not created: $NEW_SOLO"; exit 1; }

echo "✓ Outputs present"

echo "=== Comparing to baseline ==="
BAM_PASS=false; TAGS_PASS=false; SOLO_PASS=false

if compare_bams_counts "$BASELINE_BAM" "$NEW_BAM" "Baseline (sorted)" "New (unsorted+tags)"; then BAM_PASS=true; fi

echo "Checking CB/UB tags in both BAMs (first 10 records)"
BASELINE_CB_UB=$(samtools view "$BASELINE_BAM" | head -10 | grep -E "(CB:|UB:)" | wc -l || true)
NEW_CB_UB=$(samtools view "$NEW_BAM" | head -10 | grep -E "(CB:|UB:)" | wc -l || true)
echo "  Baseline CB/UB count in first 10: $BASELINE_CB_UB"
echo "  New CB/UB count in first 10:      $NEW_CB_UB"
if [[ "${BASELINE_CB_UB:-0}" -gt 0 && "${NEW_CB_UB:-0}" -gt 0 ]]; then TAGS_PASS=true; else echo "  ✗ CB/UB tags not detected in one or both"; fi

if compare_solo "$BASELINE_SOLO" "$NEW_SOLO"; then SOLO_PASS=true; fi

echo ""
echo "=== Summary ==="
echo "BAM counts:  $([ "$BAM_PASS" == true ] && echo PASSED || echo FAILED)"
echo "CB/UB tags: $([ "$TAGS_PASS" == true ] && echo PASSED || echo FAILED)"
echo "Solo:       $([ "$SOLO_PASS" == true ] && echo PASSED || echo FAILED)"

$BAM_PASS && $TAGS_PASS && $SOLO_PASS && exit 0 || exit 1
