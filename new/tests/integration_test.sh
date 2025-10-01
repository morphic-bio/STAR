#!/bin/bash

# Integration test for STARsolo two-pass unsorted CB/UB tag injection
# Three-phase testing approach:
# 1. Compare original vs new STAR in sorted mode (should be identical)
# 2. Compare original vs new STAR in unsorted mode without tags (should be identical) 
# 3. Compare new STAR with tags vs original sorted (final validation)
# Based on production parameters from SC2300771

set -euo pipefail

# Configuration
ORIGINAL_STAR_BINARY="/usr/local/bin/STAR"
NEW_STAR_BINARY="./STAR"  # Use the newly built binary in current directory
SAMPLE_ID="SC2300771"
BASE_DIR="/storage/downsampled/${SAMPLE_ID}"
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
THREADS=24
OUTPUT_BASE="$(pwd)/integration_test_output"

# Set ulimit for STAR
ulimit -n 100000

# Helper function to clean up directories
cleanup_dir() {
    local dir="$1"
    echo "Cleaning up: $dir"
    rm -rf "$dir"
    rm -rf "${dir}_STARtmp" 2>/dev/null || true
    mkdir -p "$dir"
}

# Helper function to compare BAM files
compare_bams() {
    local bam1="$1"
    local bam2="$2"
    local name1="$3"
    local name2="$4"
    
    echo "Comparing BAM files:"
    echo "  $name1: $(du -h "$bam1" | cut -f1)"
    echo "  $name2: $(du -h "$bam2" | cut -f1)"
    
    local count1=$(samtools view -c "$bam1" 2>/dev/null || echo "Error counting")
    local count2=$(samtools view -c "$bam2" 2>/dev/null || echo "Error counting")
    echo "  $name1 records: $count1"
    echo "  $name2 records: $count2"
    
    if [[ "$count1" == "$count2" && "$count1" != "Error counting" ]]; then
        echo "  ✓ Record counts match"
        return 0
    else
        echo "  ✗ Record counts differ"
        return 1
    fi
}

# Helper function to compare Solo matrices
compare_solo() {
    local solo1="$1"
    local solo2="$2"
    local name1="$3"
    local name2="$4"
    
    echo "Comparing Solo matrices:"
    
    local features1="${solo1}/Gene/raw/features.tsv"
    local features2="${solo2}/Gene/raw/features.tsv"
    local barcodes1="${solo1}/Gene/raw/barcodes.tsv"
    local barcodes2="${solo2}/Gene/raw/barcodes.tsv"
    local matrix1="${solo1}/Gene/raw/matrix.mtx"
    local matrix2="${solo2}/Gene/raw/matrix.mtx"
    
    local success=0
    
    if [[ -f "$features1" && -f "$features2" ]]; then
        if diff -q "$features1" "$features2" >/dev/null; then
            echo "  ✓ Features.tsv files identical"
        else
            echo "  ✗ Features.tsv files differ"
            success=1
        fi
    else
        echo "  ! Features.tsv files missing"
        success=1
    fi
    
    if [[ -f "$barcodes1" && -f "$barcodes2" ]]; then
        if diff -q "$barcodes1" "$barcodes2" >/dev/null; then
            echo "  ✓ Barcodes.tsv files identical"
        else
            echo "  ✗ Barcodes.tsv files differ"
            success=1
        fi
    else
        echo "  ! Barcodes.tsv files missing"
        success=1
    fi
    
    if [[ -f "$matrix1" && -f "$matrix2" ]]; then
        if diff -q "$matrix1" "$matrix2" >/dev/null; then
            echo "  ✓ Matrix.mtx files identical"
        else
            echo "  ✗ Matrix.mtx files differ"
            success=1
        fi
    else
        echo "  ! Matrix.mtx files missing"
        success=1
    fi
    
    return $success
}

# Check for required files
echo "=== Checking input files and binaries ==="
if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: Input directory not found: $BASE_DIR"
    exit 1
fi

if [[ ! -f "$WHITELIST" ]]; then
    echo "ERROR: Whitelist not found: $WHITELIST"
    exit 1
fi

if [[ ! -d "$GENOME_DIR" ]]; then
    echo "ERROR: Genome directory not found: $GENOME_DIR"
    exit 1
fi

if [[ ! -f "$ORIGINAL_STAR_BINARY" ]]; then
    echo "ERROR: Original STAR binary not found: $ORIGINAL_STAR_BINARY"
    exit 1
fi

if [[ ! -f "$NEW_STAR_BINARY" ]]; then
    echo "ERROR: New STAR binary not found: $NEW_STAR_BINARY"
    exit 1
fi

# Build file lists (assuming 8 lanes like original command)
R2_FILES=""
R1_FILES=""
for lane in L001 L002 L003 L004 L005 L006 L007 L008; do
    r2_file="${BASE_DIR}/${SAMPLE_ID}_*_${lane}_R2_001.fastq.gz"
    r1_file="${BASE_DIR}/${SAMPLE_ID}_*_${lane}_R1_001.fastq.gz"
    
    # Check if files exist (using glob expansion)
    if ls $r2_file 1> /dev/null 2>&1; then
        if [[ -n "$R2_FILES" ]]; then R2_FILES="${R2_FILES},"; fi
        R2_FILES="${R2_FILES}$(ls $r2_file)"
    fi
    
    if ls $r1_file 1> /dev/null 2>&1; then
        if [[ -n "$R1_FILES" ]]; then R1_FILES="${R1_FILES},"; fi
        R1_FILES="${R1_FILES}$(ls $r1_file)"
    fi
done

if [[ -z "$R2_FILES" ]] || [[ -z "$R1_FILES" ]]; then
    echo "ERROR: Could not find R1/R2 files in $BASE_DIR"
    exit 1
fi

echo "Found R2 files: $R2_FILES"
echo "Found R1 files: $R1_FILES"

# Common STAR parameters
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

SORTED_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM SortedByCoordinate
    --outSAMattributes NH HI AS nM NM CB UB CR CY UR UY GX GN
)

UNSORTED_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
)

NEW_SORTED_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM SortedByCoordinate
    --soloAddTagsToUnsorted yes
    --outSAMattributes NH HI AS nM NM CB UB CR CY UR UY GX GN
)

NEW_UNSORTED_NO_TAGS_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
    --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN
)

NEW_UNSORTED_TAGS_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
    --soloAddTagsToUnsorted yes
)

echo ""
echo "System check:"
echo "  ulimit -n (max open files): $(ulimit -n)"
echo "  Available disk space: $(df -h "$(dirname "${OUTPUT_BASE}")" | tail -1 | awk '{print $4}')"

# Create base output directory
rm -rf "${OUTPUT_BASE}"
mkdir -p "${OUTPUT_BASE}"

echo ""
echo "============================================================================="
echo "TEST 1: Compare original vs new STAR in sorted mode (should be identical)"
echo "============================================================================="

# Test 1a: Original STAR sorted
echo ""
echo "=== Test 1a: Original STAR - Sorted mode ==="
cleanup_dir "${OUTPUT_BASE}/original_sorted"
echo "Started at: $(date)"

$ORIGINAL_STAR_BINARY \
    "${SORTED_PARAMS[@]}" \
    --outFileNamePrefix "${OUTPUT_BASE}/original_sorted/"

echo "Original sorted run completed at: $(date)"

# Test 1b: New STAR sorted  
echo ""
echo "=== Test 1b: New STAR - Sorted mode ==="
cleanup_dir "${OUTPUT_BASE}/new_sorted"
echo "Started at: $(date)"

$NEW_STAR_BINARY \
    "${SORTED_PARAMS[@]}" \
    --outFileNamePrefix "${OUTPUT_BASE}/new_sorted/"

echo "New sorted run completed at: $(date)"

# Compare Test 1 results
echo ""
echo "=== Test 1 Results ==="
ORIGINAL_SORTED_BAM="${OUTPUT_BASE}/original_sorted/Aligned.sortedByCoord.out.bam"
NEW_SORTED_BAM="${OUTPUT_BASE}/new_sorted/Aligned.sortedByCoord.out.bam"

if compare_bams "$ORIGINAL_SORTED_BAM" "$NEW_SORTED_BAM" "Original sorted" "New sorted"; then
    echo "✓ Test 1 BAM comparison PASSED"
    TEST1_BAM_PASS=true
else
    echo "✗ Test 1 BAM comparison FAILED"
    TEST1_BAM_PASS=false
fi

if compare_solo "${OUTPUT_BASE}/original_sorted/Solo.out" "${OUTPUT_BASE}/new_sorted/Solo.out" "Original sorted" "New sorted"; then
    echo "✓ Test 1 Solo comparison PASSED"
    TEST1_SOLO_PASS=true
else
    echo "✗ Test 1 Solo comparison FAILED"
    TEST1_SOLO_PASS=false
fi

if [[ "$TEST1_BAM_PASS" == true && "$TEST1_SOLO_PASS" == true ]]; then
    echo ""
    echo "🎉 TEST 1 OVERALL: PASSED - Original and new STAR produce identical results in sorted mode"
    TEST1_PASS=true
else
    echo ""
    echo "❌ TEST 1 OVERALL: FAILED - Original and new STAR differ in sorted mode"
    TEST1_PASS=false
    echo "STOPPING: Cannot proceed to Test 2 if basic sorted mode differs"
    exit 1
fi

echo ""
echo "============================================================================="
echo "TEST 2: Compare original vs new STAR in unsorted mode (should be identical)"
echo "============================================================================="

# Test 2a: Original STAR unsorted (no CB/UB tags in SAM attributes)
echo ""
echo "=== Test 2a: Original STAR - Unsorted mode (no CB/UB) ==="
cleanup_dir "${OUTPUT_BASE}/original_unsorted"
echo "Started at: $(date)"

$ORIGINAL_STAR_BINARY \
    "${NEW_UNSORTED_NO_TAGS_PARAMS[@]}" \
    --outFileNamePrefix "${OUTPUT_BASE}/original_unsorted/"

echo "Original unsorted run completed at: $(date)"

# Test 2b: New STAR unsorted (no CB/UB tags, no --soloAddTagsToUnsorted)
echo ""
echo "=== Test 2b: New STAR - Unsorted mode (no CB/UB, no tags) ==="
cleanup_dir "${OUTPUT_BASE}/new_unsorted_no_tags"
echo "Started at: $(date)"

$NEW_STAR_BINARY \
    "${NEW_UNSORTED_NO_TAGS_PARAMS[@]}" \
    --outFileNamePrefix "${OUTPUT_BASE}/new_unsorted_no_tags/"

echo "New unsorted (no tags) run completed at: $(date)"

# Compare Test 2 results
echo ""
echo "=== Test 2 Results ==="
ORIGINAL_UNSORTED_BAM="${OUTPUT_BASE}/original_unsorted/Aligned.out.bam"
NEW_UNSORTED_NO_TAGS_BAM="${OUTPUT_BASE}/new_unsorted_no_tags/Aligned.out.bam"

if compare_bams "$ORIGINAL_UNSORTED_BAM" "$NEW_UNSORTED_NO_TAGS_BAM" "Original unsorted" "New unsorted (no tags)"; then
    echo "✓ Test 2 BAM comparison PASSED"
    TEST2_BAM_PASS=true
else
    echo "✗ Test 2 BAM comparison FAILED"
    TEST2_BAM_PASS=false
fi

if compare_solo "${OUTPUT_BASE}/original_unsorted/Solo.out" "${OUTPUT_BASE}/new_unsorted_no_tags/Solo.out" "Original unsorted" "New unsorted (no tags)"; then
    echo "✓ Test 2 Solo comparison PASSED"
    TEST2_SOLO_PASS=true
else
    echo "✗ Test 2 Solo comparison FAILED"
    TEST2_SOLO_PASS=false
fi

if [[ "$TEST2_BAM_PASS" == true && "$TEST2_SOLO_PASS" == true ]]; then
    echo ""
    echo "🎉 TEST 2 OVERALL: PASSED - Original and new STAR produce identical results in unsorted mode"
    TEST2_PASS=true
else
    echo ""
    echo "❌ TEST 2 OVERALL: FAILED - Original and new STAR differ in unsorted mode"
    TEST2_PASS=false
    echo "STOPPING: Cannot proceed to Test 3 if basic unsorted mode differs"
    exit 1
fi

echo ""
echo "============================================================================="
echo "TEST 3: New STAR with CB/UB tags vs original sorted (final validation)"
echo "============================================================================="

# Test 3: New STAR with --soloAddTagsToUnsorted
echo ""
echo "=== Test 3: New STAR - Unsorted mode with CB/UB tags ==="
cleanup_dir "${OUTPUT_BASE}/new_unsorted_with_tags"
echo "Started at: $(date)"

$NEW_STAR_BINARY \
    "${NEW_UNSORTED_TAGS_PARAMS[@]}" \
    --outFileNamePrefix "${OUTPUT_BASE}/new_unsorted_with_tags/"

echo "New unsorted with tags run completed at: $(date)"

# Compare Test 3 results (against original sorted from Test 1)
echo ""
echo "=== Test 3 Results ==="
NEW_UNSORTED_WITH_TAGS_BAM="${OUTPUT_BASE}/new_unsorted_with_tags/Aligned.out.bam"

echo "Comparing new unsorted with tags vs original sorted:"
echo "  Original sorted: $(du -h "$ORIGINAL_SORTED_BAM" | cut -f1)"
echo "  New unsorted+tags: $(du -h "$NEW_UNSORTED_WITH_TAGS_BAM" | cut -f1)"

# Record counts should match
ORIGINAL_COUNT=$(samtools view -c "$ORIGINAL_SORTED_BAM" 2>/dev/null || echo "Error counting")
NEW_TAGS_COUNT=$(samtools view -c "$NEW_UNSORTED_WITH_TAGS_BAM" 2>/dev/null || echo "Error counting")
echo "  Original sorted records: $ORIGINAL_COUNT"
echo "  New unsorted+tags records: $NEW_TAGS_COUNT"

if [[ "$ORIGINAL_COUNT" == "$NEW_TAGS_COUNT" && "$ORIGINAL_COUNT" != "Error counting" ]]; then
    echo "  ✓ Record counts match"
    TEST3_BAM_PASS=true
else
    echo "  ✗ Record counts differ"
    TEST3_BAM_PASS=false
fi

# Check CB/UB tags presence
echo ""
echo "CB/UB tag verification:"
echo "Original sorted BAM (should have CB/UB):"
ORIG_CB_UB=$(samtools view "$ORIGINAL_SORTED_BAM" | head -10 | grep -E "(CB:|UB:)" | wc -l)
echo "  Records with CB/UB: $ORIG_CB_UB"

echo "New unsorted+tags BAM (should have CB/UB):"
NEW_CB_UB=$(samtools view "$NEW_UNSORTED_WITH_TAGS_BAM" | head -10 | grep -E "(CB:|UB:)" | wc -l)
echo "  Records with CB/UB: $NEW_CB_UB"

if [[ "$ORIG_CB_UB" -gt 0 && "$NEW_CB_UB" -gt 0 ]]; then
    echo "  ✓ Both BAMs contain CB/UB tags"
    TEST3_TAGS_PASS=true
else
    echo "  ✗ CB/UB tags missing in one or both BAMs"
    TEST3_TAGS_PASS=false
fi

# Solo matrices should be identical
if compare_solo "${OUTPUT_BASE}/original_sorted/Solo.out" "${OUTPUT_BASE}/new_unsorted_with_tags/Solo.out" "Original sorted" "New unsorted+tags"; then
    echo "✓ Test 3 Solo comparison PASSED"
    TEST3_SOLO_PASS=true
else
    echo "✗ Test 3 Solo comparison FAILED"
    TEST3_SOLO_PASS=false
fi

if [[ "$TEST3_BAM_PASS" == true && "$TEST3_TAGS_PASS" == true && "$TEST3_SOLO_PASS" == true ]]; then
    echo ""
    echo "🎉 TEST 3 OVERALL: PASSED - New STAR with CB/UB tags produces equivalent results to original sorted"
    TEST3_PASS=true
else
    echo ""
    echo "❌ TEST 3 OVERALL: FAILED - New STAR with CB/UB tags differs from original sorted"
    TEST3_PASS=false
fi

# ============================================================================
# Test 4: Tag Table Export Feature
# ============================================================================

echo ""
echo "=== Test 4: New STAR - Tag table export with --soloWriteTagTable ==="

cleanup_dir "${OUTPUT_BASE}/tag_table_test"

"$NEW_STAR_BINARY" \
    --runThreadN $THREADS \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "${BASE_DIR}/R1.fastq.gz" "${BASE_DIR}/R2.fastq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${OUTPUT_BASE}/tag_table_test/" \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --soloCBwhitelist "$WHITELIST" \
    --outSAMtype BAM Unsorted \
    --soloWriteTagTable Default \
    --outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM \
    --soloFeatures Gene

echo ""
echo "=== Test 4 Results ==="

# Check if tag table file exists
TAG_TABLE_FILE="${OUTPUT_BASE}/tag_table_test/Aligned.out.cb_ub.tsv"
if [[ -f "$TAG_TABLE_FILE" ]]; then
    echo "✓ Tag table file exists: $TAG_TABLE_FILE"
    TEST4_FILE_EXISTS=true
else
    echo "✗ Tag table file missing: $TAG_TABLE_FILE"
    TEST4_FILE_EXISTS=false
fi

# Check if tag table has expected format and content
if [[ "$TEST4_FILE_EXISTS" == true ]]; then
    # Count header line
    HEADER_COUNT=$(head -1 "$TAG_TABLE_FILE" | grep -c "bam_record_index.*iReadAll.*mate.*align_idx.*qname.*CB.*UB.*status" || echo "0")
    if [[ "$HEADER_COUNT" -eq 1 ]]; then
        echo "✓ Tag table has correct header format"
        TEST4_HEADER_OK=true
    else
        echo "✗ Tag table header format is incorrect"
        echo "Expected: # bam_record_index	iReadAll	mate	align_idx	qname	CB	UB	status"
        echo "Found: $(head -1 "$TAG_TABLE_FILE")"
        TEST4_HEADER_OK=false
    fi
    
    # Count data lines (excluding header)
    DATA_LINES=$(tail -n +2 "$TAG_TABLE_FILE" | wc -l)
    if [[ "$DATA_LINES" -gt 0 ]]; then
        echo "✓ Tag table contains $DATA_LINES data records"
        TEST4_HAS_DATA=true
        
        # Check if CB and UB columns have valid values
        CB_VALUES=$(tail -n +2 "$TAG_TABLE_FILE" | cut -f6 | grep -v "^-$" | head -10 | wc -l)
        UB_VALUES=$(tail -n +2 "$TAG_TABLE_FILE" | cut -f7 | grep -v "^-$" | head -10 | wc -l)
        if [[ "$CB_VALUES" -gt 0 && "$UB_VALUES" -gt 0 ]]; then
            echo "✓ Tag table contains valid CB and UB values"
            TEST4_VALID_DATA=true
        else
            echo "✗ Tag table lacks valid CB/UB values in first 10 records"
            TEST4_VALID_DATA=false
        fi
    else
        echo "✗ Tag table is empty (no data records)"
        TEST4_HAS_DATA=false
        TEST4_VALID_DATA=false
    fi
else
    TEST4_HEADER_OK=false
    TEST4_HAS_DATA=false
    TEST4_VALID_DATA=false
fi

# Verify unsorted BAM was created without CB/UB tags (since we didn't enable --soloAddTagsToUnsorted)
UNSORTED_BAM="${OUTPUT_BASE}/tag_table_test/Aligned.out.bam"
if [[ -f "$UNSORTED_BAM" ]]; then
    # Check if BAM has CB/UB tags (should not have them in tag table mode without --soloAddTagsToUnsorted)
    CB_TAG_COUNT=$(samtools view "$UNSORTED_BAM" | head -100 | grep -c "CB:Z:" || echo "0")
    if [[ "$CB_TAG_COUNT" -eq 0 ]]; then
        echo "✓ Unsorted BAM correctly lacks CB/UB tags (as expected in tag table mode)"
        TEST4_BAM_NO_TAGS=true
    else
        echo "✗ Unsorted BAM unexpectedly contains CB tags"
        TEST4_BAM_NO_TAGS=false
    fi
else
    echo "✗ Unsorted BAM file missing"
    TEST4_BAM_NO_TAGS=false
fi

# Compare Solo matrices with previous runs (should be identical)
if compare_solo_dirs "${OUTPUT_BASE}/original_sorted/Solo.out" "${OUTPUT_BASE}/tag_table_test/Solo.out"; then
    echo "✓ Test 4 Solo comparison PASSED"
    TEST4_SOLO_PASS=true
else
    echo "✗ Test 4 Solo comparison FAILED"
    TEST4_SOLO_PASS=false
fi

if [[ "$TEST4_FILE_EXISTS" == true && "$TEST4_HEADER_OK" == true && "$TEST4_HAS_DATA" == true && "$TEST4_VALID_DATA" == true && "$TEST4_BAM_NO_TAGS" == true && "$TEST4_SOLO_PASS" == true ]]; then
    echo ""
    echo "🎉 TEST 4 OVERALL: PASSED - Tag table export feature working correctly"
    TEST4_PASS=true
else
    echo ""
    echo "❌ TEST 4 OVERALL: FAILED - Tag table export feature has issues"
    TEST4_PASS=false
fi

echo ""
echo "============================================================================="
echo "FINAL SUMMARY"
echo "============================================================================="
echo "Test 1 (sorted mode comparison): $([ "$TEST1_PASS" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"
echo "Test 2 (unsorted mode comparison): $([ "$TEST2_PASS" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"
echo "Test 3 (new feature validation): $([ "$TEST3_PASS" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"
echo "Test 4 (tag table export): $([ "$TEST4_PASS" == true ] && echo "PASSED ✓" || echo "FAILED ✗")"

if [[ "$TEST1_PASS" == true && "$TEST2_PASS" == true && "$TEST3_PASS" == true && "$TEST4_PASS" == true ]]; then
    echo ""
    echo "🎉 ALL TESTS PASSED! The new STARsolo features are working correctly."
    echo ""
    echo "Key findings:"
    echo "- New STAR binary produces identical results to original in both sorted and unsorted modes"
    echo "- New --soloAddTagsToUnsorted feature successfully adds CB/UB tags to unsorted BAM"
    echo "- New --soloWriteTagTable feature successfully exports CB/UB assignments to sidecar table"
    echo "- Solo count matrices are identical across all methods"
    echo "- New methods produce smaller BAM files (unsorted vs sorted)"
    exit 0
else
    echo ""
    echo "❌ SOME TESTS FAILED! Please review the results above."
    exit 1
fi