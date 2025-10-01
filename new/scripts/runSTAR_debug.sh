#!/bin/bash

# Debug version of runSTAR.sh with instrumentation enabled
# Tracks down memory errors and data corruption issues
# Based on debug_plan.txt implementation

set -euo pipefail

# Configuration
# Resolve absolute path to ensure the intended binary is used
NEW_STAR_BINARY="$(cd "$(dirname "$0")" && pwd)/source/STAR"                  # Use source/STAR for debug build
SAMPLE_ID="SC2300771"
BASE_DIR="/storage/JAX_sequences/${SAMPLE_ID}"
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
THREADS=8                                       # Reduced threads for debug stability
OUTPUT_BASE="/storage/Alignments"
DEBUG_LOG_DIR="/tmp/star_debug_logs"

NEW_DIR="${OUTPUT_BASE}/${SAMPLE_ID}_debug"
TEMP_DIR="/storage/tmp/${SAMPLE_ID}_debug"

# Debug-specific configuration
export STAR_DEBUG_TAG=1                        # Enable debug instrumentation
export MALLOC_CHECK_=2                         # Enable glibc malloc debugging
export MALLOC_PERTURB_=42                      # Perturb freed memory to catch use-after-free

# Args
PRE_FLIGHT_ONLY=0
for arg in "$@"; do
    case "$arg" in
        --preflight-only|--no-run)
            PRE_FLIGHT_ONLY=1
            ;;
    esac
done

# Create debug log directory
mkdir -p "$DEBUG_LOG_DIR"
DEBUG_LOG_FILE="${DEBUG_LOG_DIR}/star_debug_$(date +%Y%m%d_%H%M%S).log"

echo "=== STAR Debug Run Configuration ==="
echo "Debug binary: $NEW_STAR_BINARY"
echo "Binary SHA256: $(sha256sum "$NEW_STAR_BINARY" | awk '{print $1}')"
echo "Binary mtime:  $(stat -c %y "$NEW_STAR_BINARY")"
echo "Binary metadata strings:"
strings "$NEW_STAR_BINARY" | grep -E 'COMPILATION_TIME_PLACE|GIT_BRANCH_COMMIT_DIFF' | sed 's/^/  /' || true
echo "Sample ID: $SAMPLE_ID"
echo "Output directory: $NEW_DIR"
echo "Temp directory: $TEMP_DIR"
echo "Debug log file: $DEBUG_LOG_FILE"
echo "Threads: $THREADS (reduced for debug stability)"
echo "Debug environment:"
echo "  STAR_DEBUG_TAG=$STAR_DEBUG_TAG"
echo "  MALLOC_CHECK_=$MALLOC_CHECK_"
echo "  MALLOC_PERTURB_=$MALLOC_PERTURB_"
echo ""

# Clean up temp directories
rm -rf "$TEMP_DIR"

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

# Verify debug instrumentation is compiled in
echo "=== Verifying debug instrumentation ==="
echo "Testing debug instrumentation with quick run..."
if timeout 10 bash -c 'STAR_DEBUG_TAG=1 echo "test" | '"$NEW_STAR_BINARY"' --help 2>&1 | grep -q "DEBUG_TAG"' 2>/dev/null; then
    echo "✓ Debug instrumentation confirmed working"
elif strings "$NEW_STAR_BINARY" | grep -q "DEBUG_TAG"; then
    echo "✓ Debug instrumentation detected in binary"
else
    echo "⚠ WARNING: Debug instrumentation not detected"
    echo "  Testing with a quick run to verify..."
    # Try a more comprehensive test
    if STAR_DEBUG_TAG=1 timeout 5 "$NEW_STAR_BINARY" --version 2>&1 | grep -q "DEBUG_TAG"; then
        echo "✓ Debug instrumentation confirmed working via test run"
    else
        echo "  Consider rebuilding with: cd source && make STAR"
        echo "  Continuing anyway - instrumentation may still work..."
    fi
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

# Count files found
R2_COUNT=$(echo "$R2_FILES" | tr ',' '\n' | wc -l)
R1_COUNT=$(echo "$R1_FILES" | tr ',' '\n' | wc -l)

echo "Found $R2_COUNT R2 files and $R1_COUNT R1 files"
echo "R2: $R2_FILES"
echo "R1: $R1_FILES"

if [[ $R2_COUNT -eq 1 && $R1_COUNT -eq 1 ]]; then
    echo "ℹ Single lane detected - this is fine for debugging"
elif [[ $R2_COUNT -ne $R1_COUNT ]]; then
    echo "⚠ WARNING: Mismatch between R2 ($R2_COUNT) and R1 ($R1_COUNT) file counts"
fi

# Common STAR params (debug-optimized)
COMMON_PARAMS=(
    --runThreadN $THREADS
    --outTmpDir "$TEMP_DIR"
    --soloType CB_UMI_Simple
    --soloCBlen 16
    --soloUMIlen 12
    --soloUMIstart 17
    --soloCBstart 1
    --soloBarcodeReadLength 0
    --soloCBwhitelist "$WHITELIST"
    --genomeDir "$GENOME_DIR"
    --limitIObufferSize 50000000 50000000        # Reduced buffer sizes for debug
    --outSJtype None
    --outBAMcompression 6                        # Reduced compression for faster debug
    --soloMultiMappers Unique
    --alignIntronMax 1
    --alignMatesGapMax 0
    --outFilterMismatchNmax 6
    --outFilterMismatchNoverReadLmax 1.0
    --outFilterMatchNmin 25
    --outSAMunmapped None
    --outFilterMatchNminOverLread 0
    --outFilterMultimapNmax 10000
    --outFilterMultimapScoreRange 1
    --outSAMmultNmax 10000
    --winAnchorMultimapNmax 200
    --outSAMprimaryFlag AllBestScore
    --outFilterScoreMin 0
    --outFilterScoreMinOverLread 0
    --outSAMattributes NH HI AS nM NM CB UB CR CY UR UY GX GN  # Added CB UB for debug
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

NEW_TAG_TABLE_PARAMS=(
    "${COMMON_PARAMS[@]}"
    --outSAMtype BAM Unsorted
    --soloWriteTagTable Default
    --soloAddTagsToUnsorted yes                  # Enable both modes for comprehensive debug
)

# Function to monitor system resources
monitor_resources() {
    local pid=$1
    local log_file=$2
    echo "=== Resource monitoring started for PID $pid ===" >> "$log_file"
    while kill -0 "$pid" 2>/dev/null; do
        echo "$(date): Memory usage:" >> "$log_file"
        ps -p "$pid" -o pid,ppid,vsz,rss,pmem,pcpu,time,cmd >> "$log_file" 2>/dev/null || break
        echo "$(date): System memory:" >> "$log_file"
        free -h >> "$log_file"
        echo "---" >> "$log_file"
        sleep 30
    done
    echo "=== Resource monitoring ended ===" >> "$log_file"
}

# Function to analyze debug output in real-time
analyze_debug_output() {
    local log_file=$1
    echo "=== Debug Output Analysis ==="
    
    # Count debug messages by type
    local debug_count
    local error_count
    local warning_count
    debug_count=$(grep -c "DEBUG_TAG" "$log_file" 2>/dev/null || true)
    error_count=$(grep -c "ERROR:" "$log_file" 2>/dev/null || true)
    warning_count=$(grep -c "WARNING:" "$log_file" 2>/dev/null || true)
    [[ -n "$debug_count" ]] || debug_count=0
    [[ -n "$error_count" ]] || error_count=0
    [[ -n "$warning_count" ]] || warning_count=0
    
    echo "Debug messages: $debug_count"
    echo "Errors: $error_count"
    echo "Warnings: $warning_count"
    
    if [[ $error_count -gt 0 ]]; then
        echo ""
        echo "=== ERRORS DETECTED ==="
        grep "ERROR:" "$log_file" | tail -10
    fi
    
    if [[ $warning_count -gt 0 ]]; then
        echo ""
        echo "=== WARNINGS DETECTED ==="
        grep "WARNING:" "$log_file" | tail -5
    fi
    
    # Look for specific memory issues
    if grep -q "Invalid read index" "$log_file" 2>/dev/null; then
        echo ""
        echo "=== MEMORY CORRUPTION DETECTED ==="
        grep -A 3 -B 3 "Invalid read index" "$log_file"
        if ! grep -q "trailer: full=" "$log_file" 2>/dev/null; then
            echo ""
            echo "Hint: Missing 'trailer: full=' debug lines suggests an outdated binary was used."
            echo "      Rebuild and rerun: (cd source && make STAR)"
        fi
    fi
    
    if grep -q "Allocation failed" "$log_file" 2>/dev/null; then
        echo ""
        echo "=== ALLOCATION FAILURE DETECTED ==="
        grep -A 3 -B 3 "Allocation failed" "$log_file"
    fi
}

# Pre-flight checks
echo "=== Pre-flight system checks ==="
echo "Available memory:"
free -h
echo "Available disk space:"
df -h "$OUTPUT_BASE"
# Create temp directory before checking disk space
mkdir -p "$TEMP_DIR"
df -h "$TEMP_DIR"
echo "Open file limit: $(ulimit -n)"
echo ""
rm -rf "$TEMP_DIR"

if [[ $PRE_FLIGHT_ONLY -eq 1 ]]; then
    echo "Preflight-only mode requested; exiting before STAR run."
    exit 0
fi

# Run new STAR with debug instrumentation
cleanup_dir "$NEW_DIR"

echo "=== Running STAR with debug instrumentation enabled ==="
echo "Command line:"
echo "$NEW_STAR_BINARY \\"
printf '    %s \\\n' "${NEW_TAG_TABLE_PARAMS[@]}"
echo "    --outFileNamePrefix \"${NEW_DIR}/\""
echo ""
echo "Starting STAR run at $(date)..."
echo "Debug output will be logged to: $DEBUG_LOG_FILE"
echo "Monitor with: tail -f $DEBUG_LOG_FILE"
echo ""

# Start STAR with debug logging
{
    echo "=== STAR Debug Run Started at $(date) ==="
    echo "Command: $NEW_STAR_BINARY ${NEW_TAG_TABLE_PARAMS[*]} --outFileNamePrefix ${NEW_DIR}/"
    echo "Environment: STAR_DEBUG_TAG=$STAR_DEBUG_TAG MALLOC_CHECK_=$MALLOC_CHECK_ MALLOC_PERTURB_=$MALLOC_PERTURB_"
    echo "Working directory: $(pwd)"
    echo "=== Debug Output ==="
} > "$DEBUG_LOG_FILE"

# Run STAR and capture both stdout and stderr
set +e  # Don't exit on STAR failure - we want to analyze the logs
"$NEW_STAR_BINARY" \
    "${NEW_TAG_TABLE_PARAMS[@]}" \
    --outFileNamePrefix "${NEW_DIR}/" \
    >> "$DEBUG_LOG_FILE" 2>&1 &

STAR_PID=$!
echo "STAR process started with PID: $STAR_PID"

# Start resource monitoring in background
monitor_resources $STAR_PID "${DEBUG_LOG_FILE}.resources" &
MONITOR_PID=$!

# Wait for STAR to complete
wait $STAR_PID
STAR_EXIT_CODE=$?

# Stop resource monitoring
kill $MONITOR_PID 2>/dev/null || true
wait $MONITOR_PID 2>/dev/null || true

echo ""
echo "=== STAR Debug Run Completed at $(date) ==="
echo "Exit code: $STAR_EXIT_CODE"

# Analyze the debug output
analyze_debug_output "$DEBUG_LOG_FILE"

# Check output files
echo ""
echo "=== Output File Analysis ==="
if [[ -f "${NEW_DIR}/Aligned.out.cb_ub.bin" ]]; then
    bin_size=$(stat -c%s "${NEW_DIR}/Aligned.out.cb_ub.bin")
    echo "✓ Binary tag table created: ${NEW_DIR}/Aligned.out.cb_ub.bin ($bin_size bytes)"
else
    echo "✗ Binary tag table NOT created"
fi

if [[ -f "${NEW_DIR}/Aligned.out.bam" ]]; then
    bam_size=$(stat -c%s "${NEW_DIR}/Aligned.out.bam")
    echo "✓ BAM file created: ${NEW_DIR}/Aligned.out.bam ($bam_size bytes)"
else
    echo "✗ BAM file NOT created"
fi

# Summary
echo ""
echo "=== Debug Run Summary ==="
echo "STAR exit code: $STAR_EXIT_CODE"
echo "Debug log: $DEBUG_LOG_FILE"
echo "Resource log: ${DEBUG_LOG_FILE}.resources"
echo "Output directory: $NEW_DIR"

if [[ $STAR_EXIT_CODE -eq 0 ]]; then
    echo "✅ STAR completed successfully"
else
    echo "❌ STAR failed with exit code $STAR_EXIT_CODE"
    echo ""
    echo "=== Last 20 lines of debug log ==="
    tail -20 "$DEBUG_LOG_FILE"
fi

echo ""
echo "To analyze the full debug log:"
echo "  less $DEBUG_LOG_FILE"
echo "To search for specific issues:"
echo "  grep -i 'error\\|warning\\|invalid\\|allocation' $DEBUG_LOG_FILE"
echo "To view resource usage:"
echo "  less ${DEBUG_LOG_FILE}.resources"

exit $STAR_EXIT_CODE
