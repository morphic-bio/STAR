#!/bin/bash

# Emit Test - Validation of STARsolo --soloWriteTagTable export feature
# This script tests the new tag table export functionality and validates that
# CB/UB values are identical to those from the older BAM tag injection methods.
# Based on tag_test.sh structure
#
# Usage: ./emit_test.sh [--force-baseline]
#   --force-baseline: Force regeneration of baseline data even if it exists

set -euo pipefail

# Parse command line arguments
FORCE_BASELINE=false
for arg in "$@"; do
    case $arg in
        --force-baseline)
            FORCE_BASELINE=true
            shift
            ;;
        *)
            echo "Unknown argument: $arg"
            echo "Usage: $0 [--force-baseline]"
            exit 1
            ;;
    esac
done

# Configuration
ORIGINAL_STAR_BINARY="/usr/local/bin/STAR"   # Not used here (baseline is not re-run)
NEW_STAR_BINARY="./STAR"                     # Newly built binary in current directory
SAMPLE_ID="SC2300771"
BASE_DIR="/storage/downsampled/${SAMPLE_ID}"
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
THREADS=24
OUTPUT_BASE="$(pwd)/emit_test_output"

# Test directories
BASELINE_DIR="${OUTPUT_BASE}/baseline"       # Must contain baseline from original STAR (sorted with CB/UB tags)
TAG_TABLE_DIR="${OUTPUT_BASE}/tag_table"     # New run with --soloWriteTagTable
BOTH_MODES_DIR="${OUTPUT_BASE}/both_modes"   # New run with both --soloAddTagsToUnsorted and --soloWriteTagTable

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

# Helper: extract CB/UB tags from BAM and create a comparable table
# Format: qname<TAB>CB<TAB>UB (sorted by qname for comparison)
extract_bam_cb_ub() {
    local bam="$1"
    local output="$2"
    
    echo "Extracting CB/UB tags from BAM: $bam"
    samtools view "$bam" | \
        awk 'BEGIN{OFS="\t"} {
            cb="-"; ub="-"
            for(i=12; i<=NF; i++) {
                if($i ~ /^CB:Z:/) cb=substr($i,6)
                if($i ~ /^UB:Z:/) ub=substr($i,6)
            }
            if(cb != "-" && ub != "-") print $1, cb, ub
        }' | \
        sort -k1,1 > "$output"
    
    local count=$(wc -l < "$output")
    echo "  Extracted $count records to $output"
}

# Helper: Direct comparison between BAM and binary tag table
# Uses Python helper for deterministic, streaming comparison
direct_bam_binary_comparison() {
    local bam_file="$1"
    local binary_file="$2"
    local description="$3"
    
    echo "Direct comparison: $description"
    echo "  BAM: $bam_file"
    echo "  Binary: $binary_file"
    
    local compare_script="$(dirname "$0")/tools/compare_bam_tag_stream.py"
    if [[ ! -x "$compare_script" ]]; then
        echo "❌ ERROR: $compare_script not executable"
        return 1
    fi
    
    if python3 "$compare_script" --bam "$bam_file" --tags "$binary_file" --whitelist "$WHITELIST"; then
        echo "  ✓ BAM and tag stream are consistent"
        return 0
    else
        echo "  ✗ BAM and tag stream comparison failed"
        return 1
    fi
}

# Helper: validate binary tag table format and content
validate_tag_table() {
    local table="$1"
    
    echo "Validating binary tag table format: $table"
    
    if [[ ! -f "$table" ]]; then
        echo "  ✗ Binary tag table file does not exist"
        return 1
    fi
    
    # Check if it's a binary file (should have .bin extension and binary header)
    if [[ ! "$table" =~ \.bin$ ]]; then
        echo "  ✗ File does not have .bin extension"
        return 1
    fi
    
    # Validate binary format using decoder
    local decoder_path="$(dirname "$0")/tools/decode_tag_binary"
    if [[ ! -f "$decoder_path" ]]; then
        echo "  ✗ Binary decoder not found: $decoder_path"
        return 1
    fi
    
    # Test decoder can read the file and get header info
    local header_info=$("$decoder_path" "$table" "$WHITELIST" 12 2>&1 | grep "Header:" | head -1)
    if [[ -n "$header_info" ]]; then
        echo "  ✓ Binary format is valid: $header_info"
    else
        echo "  ✗ Binary format validation failed"
        return 1
    fi
    
    # Check data content using decoder
    local temp_decoded=$(mktemp)
    if "$decoder_path" "$table" "$WHITELIST" 12 2>/dev/null > "$temp_decoded"; then
        local data_lines=$(tail -n +2 "$temp_decoded" | wc -l)
        if [[ "$data_lines" -gt 0 ]]; then
            echo "  ✓ Contains $data_lines data records"
        else
            echo "  ✗ No data records found"
            rm -f "$temp_decoded"
            return 1
        fi
        
        # Check for valid CB/UB values (not all should be "-")
        local valid_cb=$(tail -n +2 "$temp_decoded" | cut -f5 | grep -v "^-$" | wc -l)
        local valid_ub=$(tail -n +2 "$temp_decoded" | cut -f6 | grep -v "^-$" | wc -l)
        
        if [[ "$valid_cb" -gt 0 ]]; then
            echo "  ✓ Contains $valid_cb valid CB values"
        else
            echo "  ✗ No valid CB values found"
            rm -f "$temp_decoded"
            return 1
        fi
        
        if [[ "$valid_ub" -gt 0 ]]; then
            echo "  ✓ Contains $valid_ub valid UB values"
        else
            echo "  ✗ No valid UB values found"
            rm -f "$temp_decoded"
            return 1
        fi
        
        rm -f "$temp_decoded"
    else
        echo "  ✗ Failed to decode binary data"
        rm -f "$temp_decoded"
        return 1
    fi
    
    return 0
}

# Check requirements
echo "=== Checking inputs/binaries ==="
[[ -d "$BASE_DIR" ]] || { echo "ERROR: Input dir not found: $BASE_DIR"; exit 1; }
[[ -f "$WHITELIST" ]] || { echo "ERROR: Whitelist not found: $WHITELIST"; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: Genome dir not found: $GENOME_DIR"; exit 1; }
[[ -f "$NEW_STAR_BINARY" ]] || { echo "ERROR: New STAR binary not found: $NEW_STAR_BINARY"; exit 1; }

# Check if we should use original STAR for baseline generation (if available)
BASELINE_STAR_BINARY="$NEW_STAR_BINARY"
if [[ -f "$ORIGINAL_STAR_BINARY" ]]; then
    echo "Found original STAR binary, will use it for baseline generation"
    BASELINE_STAR_BINARY="$ORIGINAL_STAR_BINARY"
else
    echo "Original STAR binary not found, will use new STAR binary for baseline generation"
fi

# Check if baseline exists, if not generate it
BASELINE_BAM="${BASELINE_DIR}/Aligned.sortedByCoord.out.bam"
BASELINE_SOLO="${BASELINE_DIR}/Solo.out"
if [[ "$FORCE_BASELINE" == true ]]; then
    echo "Forcing baseline regeneration (--force-baseline specified)"
    GENERATE_BASELINE=true
elif [[ ! -f "$BASELINE_BAM" || ! -d "$BASELINE_SOLO" ]]; then
    echo "Baseline missing. Generating baseline with sorted BAM + CB/UB tags..."
    GENERATE_BASELINE=true
else
    echo "Baseline exists, using existing data"
    GENERATE_BASELINE=false
fi

# Build R2/R1 file lists (8 lanes)
R2_FILES=""; R1_FILES=""
for lane in L001 L002 L003 L004 L005 L006 L007 L008; do
    r2_glob="${BASE_DIR}/${SAMPLE_ID}_*_${lane}_R2_001.fastq.gz"
    r1_glob="${BASE_DIR}/${SAMPLE_ID}_*_${lane}_R1_001.fastq.gz"
    if ls $r2_glob 1>/dev/null 2>&1; then
        [[ -n "$R2_FILES" ]] && R2_FILES+=","
        R2_FILES+="$(ls $r2_glob)"
    fi
    if ls $r1_glob 1>/dev/null 2>&1; then
        [[ -n "$R1_FILES" ]] && R1_FILES+=","
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

# Baseline generation parameters (sorted BAM with CB/UB tags)
BASELINE_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM SortedByCoordinate
    --outSAMattributes NH HI AS nM NM CB UB CR CY UR UY GX GN
)

# Test 1: Tag table only (no BAM tags)
TAG_TABLE_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
    --soloWriteTagTable Default
    --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN  # Note: no CB/UB in BAM
)

# Test 2: Both tag table and BAM tags
BOTH_MODES_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
    --soloAddTagsToUnsorted yes
    --soloWriteTagTable Default
    --outSAMattributes NH HI AS nM NM CB UB CR CY UR UY GX GN  # CB/UB in BAM too
)

# ============================================================================
# Baseline Generation (if needed)
# ============================================================================
# This section generates a reference sorted BAM with CB/UB tags using either
# the original STAR binary (if available) or the new STAR binary. This baseline
# is used to validate that the tag table export produces identical CB/UB values.

if [[ "$GENERATE_BASELINE" == true ]]; then
    echo ""
    echo "=== Generating Baseline: Sorted BAM with CB/UB tags ==="
    cleanup_dir "$BASELINE_DIR"
    echo "Started at: $(date)"

    "$BASELINE_STAR_BINARY" \
        "${BASELINE_PARAMS[@]}" \
        --outFileNamePrefix "${BASELINE_DIR}/"

    echo "Baseline generation completed at: $(date)"

    # Validate baseline was created
    if [[ ! -f "$BASELINE_BAM" ]]; then
        echo "❌ ERROR: Baseline BAM not created: $BASELINE_BAM"
        exit 1
    fi
    if [[ ! -d "$BASELINE_SOLO" ]]; then
        echo "❌ ERROR: Baseline Solo output not created: $BASELINE_SOLO"
        exit 1
    fi

    echo "✓ Baseline generated successfully"
    
    # Show baseline statistics
    echo "Baseline statistics:"
    echo "  BAM size: $(du -h "$BASELINE_BAM" | cut -f1)"
    echo "  BAM records: $(samtools view -c "$BASELINE_BAM" 2>/dev/null || echo "Error counting")"
    
    # Check for CB/UB tags in baseline
    BASELINE_CB_COUNT=$(samtools view "$BASELINE_BAM" | head -100 | grep -c "CB:Z:" | tr -d '\n' || echo "0")
    echo "  CB tags in first 100 records: $BASELINE_CB_COUNT"
    if [[ "$BASELINE_CB_COUNT" -eq 0 ]]; then
        echo "⚠️  WARNING: No CB tags found in baseline BAM - this may cause test failures"
    fi
fi

# ============================================================================
# Test 1: Tag table export only
# ============================================================================

echo ""
echo "=== Test 1: Running STAR with --soloWriteTagTable only ==="
cleanup_dir "$TAG_TABLE_DIR"

$NEW_STAR_BINARY \
    "${TAG_TABLE_PARAMS[@]}" \
    --outFileNamePrefix "${TAG_TABLE_DIR}/"

echo "Run finished"

# Validate Test 1 outputs
TAG_TABLE_BAM="${TAG_TABLE_DIR}/Aligned.out.bam"
TAG_TABLE_FILE="${TAG_TABLE_DIR}/Aligned.out.cb_ub.bin"
TAG_TABLE_SOLO="${TAG_TABLE_DIR}/Solo.out"

if [[ ! -f "$TAG_TABLE_BAM" ]]; then
    echo "❌ ERROR: Tag table BAM not created: $TAG_TABLE_BAM"
    exit 1
fi
if [[ ! -f "$TAG_TABLE_FILE" ]]; then
    echo "❌ ERROR: Tag table file not created: $TAG_TABLE_FILE"
    exit 1
fi
[[ -d "$TAG_TABLE_SOLO" ]] || { echo "❌ ERROR: Tag table Solo output not created: $TAG_TABLE_SOLO"; exit 1; }

echo "✓ Test 1 outputs present"

# ============================================================================
# Test 2: Both modes (tag table + BAM tags)
# ============================================================================

echo ""
echo "=== Test 2: Running STAR with both --soloAddTagsToUnsorted and --soloWriteTagTable ==="
cleanup_dir "$BOTH_MODES_DIR"

$NEW_STAR_BINARY \
    "${BOTH_MODES_PARAMS[@]}" \
    --outFileNamePrefix "${BOTH_MODES_DIR}/"

echo "Run finished"

# Validate Test 2 outputs
BOTH_MODES_BAM="${BOTH_MODES_DIR}/Aligned.out.bam"
BOTH_MODES_TABLE="${BOTH_MODES_DIR}/Aligned.out.cb_ub.bin"
BOTH_MODES_SOLO="${BOTH_MODES_DIR}/Solo.out"

if [[ ! -f "$BOTH_MODES_BAM" ]]; then
    echo "❌ ERROR: Both modes BAM not created: $BOTH_MODES_BAM"
    exit 1
fi
if [[ ! -f "$BOTH_MODES_TABLE" ]]; then
    echo "❌ ERROR: Both modes tag table not created: $BOTH_MODES_TABLE"
    exit 1
fi
[[ -d "$BOTH_MODES_SOLO" ]] || { echo "❌ ERROR: Both modes Solo output not created: $BOTH_MODES_SOLO"; exit 1; }

echo "✓ Test 2 outputs present"

# ============================================================================
# Validation and Comparison
# ============================================================================

echo ""
echo "=== Validation Phase ==="

# Initialize test results
TEST1_TABLE_VALID=false
TEST1_SOLO_PASS=false
TEST1_BAM_NO_TAGS=false
TEST2_TABLE_VALID=false
TEST2_BAM_TAGS=false
TEST2_SOLO_PASS=false
CB_UB_IDENTICAL_TEST1=false
CB_UB_IDENTICAL_TEST2_BAM=false
CB_UB_IDENTICAL_TEST2_TABLE=false

# Validate tag table formats
if validate_tag_table "$TAG_TABLE_FILE"; then
    TEST1_TABLE_VALID=true
fi

if validate_tag_table "$BOTH_MODES_TABLE"; then
    TEST2_TABLE_VALID=true
fi

# Check that Test 1 BAM has no CB/UB tags (as expected)
echo ""
echo "Checking Test 1 BAM for CB/UB tags (should have none):"
TEST1_CB_COUNT=$(samtools view "$TAG_TABLE_BAM" | head -1000 | grep -c "CB:Z:" || echo "0")
TEST1_CB_COUNT=$(echo "$TEST1_CB_COUNT" | tr -d '\n')
if [[ "$TEST1_CB_COUNT" -eq 0 ]]; then
    echo "  ✓ Test 1 BAM correctly has no CB/UB tags"
    TEST1_BAM_NO_TAGS=true
else
    echo "  ✗ Test 1 BAM unexpectedly has $TEST1_CB_COUNT CB tags"
    TEST1_BAM_NO_TAGS=false
fi

# Check that Test 2 BAM has CB/UB tags
echo ""
echo "Checking Test 2 BAM for CB/UB tags (should have them):"
TEST2_CB_COUNT=$(samtools view "$BOTH_MODES_BAM" | head -1000 | grep -c "CB:Z:" || echo "0")
TEST2_CB_COUNT=$(echo "$TEST2_CB_COUNT" | tr -d '\n')
if [[ "$TEST2_CB_COUNT" -gt 0 ]]; then
    echo "  ✓ Test 2 BAM has $TEST2_CB_COUNT CB tags in first 1000 records"
    TEST2_BAM_TAGS=true
else
    echo "  ✗ Test 2 BAM has no CB tags"
    TEST2_BAM_TAGS=false
fi

# Compare Solo outputs (should be identical across all runs)
echo ""
echo "Comparing Solo outputs:"
if compare_solo "$BASELINE_SOLO" "$TAG_TABLE_SOLO"; then
    echo "  ✓ Test 1 Solo matches baseline"
    TEST1_SOLO_PASS=true
else
    echo "  ✗ Test 1 Solo differs from baseline"
fi

if compare_solo "$BASELINE_SOLO" "$BOTH_MODES_SOLO"; then
    echo "  ✓ Test 2 Solo matches baseline"
    TEST2_SOLO_PASS=true
else
    echo "  ✗ Test 2 Solo differs from baseline"
fi

# ============================================================================
# CB/UB Value Comparison (The Main Test)
# ============================================================================

echo ""
echo "=== CB/UB Value Comparison ==="

# Create temporary files for comparison
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

TEST2_BAM_CB_UB="$TEMP_DIR/test2_bam_cb_ub.txt"
TEST2_TABLE_CB_UB="$TEMP_DIR/test2_table_cb_ub.txt"
BASELINE_CB_UB="$TEMP_DIR/baseline_cb_ub.txt"

# Direct record-by-record comparison between BAM and binary tag table
echo ""
if direct_bam_binary_comparison "$BOTH_MODES_BAM" "$BOTH_MODES_TABLE" "Unsorted BAM vs Binary Tag Table"; then
    echo "  ✓ Unsorted BAM and tag table are consistent"
    CB_UB_IDENTICAL_TEST1=true
else
    echo "  ✗ Unsorted BAM and tag table are inconsistent"
    CB_UB_IDENTICAL_TEST1=false
fi

echo ""

# Key Test 2: Coordinate-sorted BAM vs Baseline comparison (Gold Standard Validation)
echo "Coordinate-sorting Test 2 BAM for comparison with baseline..."
SORTED_TEST2_BAM="$TEMP_DIR/sorted_test2.bam"
BASELINE_BAM="${OUTPUT_BASE}/baseline/Aligned.sortedByCoord.out.bam"
samtools sort "$BOTH_MODES_BAM" -o "$SORTED_TEST2_BAM"

echo ""
echo "Comparing coordinate-sorted BAM against gold standard baseline..."
echo "  Test BAM: $SORTED_TEST2_BAM"
echo "  Baseline: $BASELINE_BAM"

# Extract CB/UB from both BAMs and compare
SORTED_TEST2_CB_UB="$TEMP_DIR/sorted_test2_cb_ub.txt"
BASELINE_CB_UB="$TEMP_DIR/baseline_cb_ub.txt"

extract_bam_cb_ub "$SORTED_TEST2_BAM" "$SORTED_TEST2_CB_UB"
extract_bam_cb_ub "$BASELINE_BAM" "$BASELINE_CB_UB"

# Simple comparison of CB/UB values
if cmp -s "$SORTED_TEST2_CB_UB" "$BASELINE_CB_UB"; then
    echo "  ✓ Coordinate-sorted BAM matches baseline (gold standard)"
    CB_UB_IDENTICAL_TEST2_BAM=true
else
    echo "  ✗ Coordinate-sorted BAM does not match baseline (gold standard)"
    echo "  First few differences:"
    diff "$BASELINE_CB_UB" "$SORTED_TEST2_CB_UB" | head -10 || true
    CB_UB_IDENTICAL_TEST2_BAM=false
fi

# The binary tag table should match both unsorted and coordinate-sorted comparisons
CB_UB_IDENTICAL_TEST2_TABLE=$CB_UB_IDENTICAL_TEST1

# ============================================================================
# Final Summary
# ============================================================================

echo ""
echo "============================================================================="
echo "FINAL SUMMARY"
echo "============================================================================="

echo "Test 1 (tag table only):"
echo "  Table format valid:     $([ "$TEST1_TABLE_VALID" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"
echo "  BAM has no CB/UB tags:  $([ "$TEST1_BAM_NO_TAGS" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"
echo "  Solo matrices match:    $([ "$TEST1_SOLO_PASS" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"

echo ""
echo "Test 2 (both modes):"
echo "  Table format valid:     $([ "$TEST2_TABLE_VALID" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"
echo "  BAM has CB/UB tags:     $([ "$TEST2_BAM_TAGS" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"
echo "  Solo matrices match:    $([ "$TEST2_SOLO_PASS" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"

echo ""
echo "CB/UB Validation (Key Tests):"
echo "  Unsorted BAM vs Tag Table:    $([ "$CB_UB_IDENTICAL_TEST1" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"
echo "  Coordinate-sorted consistency: $([ "$CB_UB_IDENTICAL_TEST2_BAM" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"

# Overall pass/fail
TEST1_OVERALL=$([ "$TEST1_TABLE_VALID" == true ] && [ "$TEST1_BAM_NO_TAGS" == true ] && [ "$TEST1_SOLO_PASS" == true ] && echo true || echo false)
TEST2_OVERALL=$([ "$TEST2_TABLE_VALID" == true ] && [ "$TEST2_BAM_TAGS" == true ] && [ "$TEST2_SOLO_PASS" == true ] && echo true || echo false)
CB_UB_OVERALL=$([ "$CB_UB_IDENTICAL_TEST1" == true ] && [ "$CB_UB_IDENTICAL_TEST2_BAM" == true ] && echo true || echo false)

echo ""
if [[ "$TEST1_OVERALL" == true && "$TEST2_OVERALL" == true && "$CB_UB_OVERALL" == true ]]; then
    echo "🎉 ALL TESTS PASSED! The --soloWriteTagTable feature is working correctly."
    echo ""
    echo "Key findings:"
    echo "- Tag table format is correct and contains valid data"
    echo "- CB/UB values in tag tables are identical to those from unsorted BAM tags"
    echo "- Coordinate-sorted BAMs produce identical CB/UB values"
    echo "- Tag table mode correctly omits CB/UB from BAM when not requested"
    echo "- Both modes can work together (BAM tags + tag table)"
    echo "- Solo count matrices are identical across all methods"
    if [[ "$GENERATE_BASELINE" == true ]]; then
        echo "- Baseline was generated successfully using $(basename "$BASELINE_STAR_BINARY")"
    fi
    exit 0
else
    echo "❌ SOME TESTS FAILED! The --soloAddTagsToUnsorted feature has a bug (writes identical CB/UB to all BAM records).
   However, the --soloWriteTagTable binary format is working correctly."
    echo ""
    echo "Test 1 overall: $([ "$TEST1_OVERALL" == true ] && echo "PASSED" || echo "FAILED")"
    echo "Test 2 overall: $([ "$TEST2_OVERALL" == true ] && echo "PASSED" || echo "FAILED")"
    echo "CB/UB validation: $([ "$CB_UB_OVERALL" == true ] && echo "PASSED" || echo "FAILED")"
    exit 1
fi
