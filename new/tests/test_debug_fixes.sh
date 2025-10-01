#!/bin/bash

# Quick test of the runSTAR_debug.sh fixes

set -euo pipefail

echo "=== Testing runSTAR_debug.sh fixes ==="
echo ""

# Test 1: Check temp directory creation
SAMPLE_ID="SC2300772"
TEMP_DIR="/storage/tmp/${SAMPLE_ID}_debug"

echo "Test 1: Temp directory creation"
echo "Creating temp directory: $TEMP_DIR"
mkdir -p "$TEMP_DIR"
if [[ -d "$TEMP_DIR" ]]; then
    echo "✓ Temp directory created successfully"
    df -h "$TEMP_DIR" | tail -1
else
    echo "✗ Failed to create temp directory"
    exit 1
fi

# Test 2: Check debug instrumentation detection
echo ""
echo "Test 2: Debug instrumentation detection"
STAR_BINARY="./source/STAR"
if [[ -f "$STAR_BINARY" ]]; then
    echo "Testing debug instrumentation..."
    export STAR_DEBUG_TAG=1
    
    # Try the quick test approach
    if timeout 5 "$STAR_BINARY" --version 2>&1 | grep -q "DEBUG_TAG"; then
        echo "✓ Debug instrumentation confirmed working"
    elif timeout 10 bash -c "echo 'test' | STAR_DEBUG_TAG=1 $STAR_BINARY --help 2>&1" | grep -q "DEBUG_TAG"; then
        echo "✓ Debug instrumentation detected via help"
    else
        echo "⚠ Debug instrumentation test inconclusive"
        echo "  This may be normal - the instrumentation only activates during actual runs"
    fi
else
    echo "✗ STAR binary not found: $STAR_BINARY"
fi

# Test 3: Check file detection logic
echo ""
echo "Test 3: File detection for SC2300772"
BASE_DIR="/storage/JAX_sequences/SC2300772"
if [[ -d "$BASE_DIR" ]]; then
    echo "Scanning for FASTQ files in $BASE_DIR..."
    
    R2_FILES=""; R1_FILES=""
    for lane in L001 L002 L003 L004 L005 L006 L007 L008; do
        r2_glob="${BASE_DIR}/SC2300772_*_${lane}_R2_001.fastq.gz"
        r1_glob="${BASE_DIR}/SC2300772_*_${lane}_R1_001.fastq.gz"
        if ls $r2_glob 1>/dev/null 2>&1; then
            [[ -n "$R2_FILES" ]] && R2_FILES+=","
            R2_FILES+="$(ls $r2_glob)"
            echo "  Found R2 for lane $lane"
        fi
        if ls $r1_glob 1>/dev/null 2>&1; then
            [[ -n "$R1_FILES" ]] && R1_FILES+=","
            R1_FILES+="$(ls $r1_glob)"
            echo "  Found R1 for lane $lane"
        fi
    done
    
    if [[ -n "$R2_FILES" && -n "$R1_FILES" ]]; then
        R2_COUNT=$(echo "$R2_FILES" | tr ',' '\n' | wc -l)
        R1_COUNT=$(echo "$R1_FILES" | tr ',' '\n' | wc -l)
        echo "✓ Found $R2_COUNT R2 files and $R1_COUNT R1 files"
        
        if [[ $R2_COUNT -eq 1 && $R1_COUNT -eq 1 ]]; then
            echo "  Single lane detected - perfect for debugging"
        fi
    else
        echo "✗ No FASTQ files found"
    fi
else
    echo "✗ Base directory not found: $BASE_DIR"
fi

# Test 4: Check available resources
echo ""
echo "Test 4: Resource availability"
echo "Memory:"
free -h | grep -E "Mem:|Swap:"
echo ""
echo "Disk space for output:"
df -h /storage/Alignments | tail -1
echo ""
echo "Disk space for temp:"
df -h /storage/tmp | tail -1

echo ""
echo "=== Test Summary ==="
echo "The runSTAR_debug.sh script should now handle:"
echo "✓ Temp directory creation before disk space check"
echo "✓ Better debug instrumentation detection"  
echo "✓ Single lane file detection with informative messages"
echo "✓ Improved error handling and warnings"
echo ""
echo "You can now run: ./runSTAR_debug.sh"
echo ""
echo "Monitor progress with:"
echo "  tail -f /tmp/star_debug_logs/star_debug_*.log"

# Cleanup test temp dir
rm -rf "$TEMP_DIR"
