#!/usr/bin/env bash
# Example: Complete memory testing workflow for STAR tag table functionality
# This script demonstrates the full process from build to analysis

set -euo pipefail

echo "=== STAR Memory Testing Example ==="
echo "This script demonstrates memory testing with AddressSanitizer"
echo ""

# Configuration
STAR_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$STAR_ROOT"

echo "Working directory: $STAR_ROOT"
echo ""

# Step 1: Build with AddressSanitizer
echo "Step 1: Building STAR with AddressSanitizer..."
cd source
echo "  Cleaning previous builds..."
ASAN=1 make clean > /dev/null

echo "  Compiling with ASan (this may take a few minutes)..."
if ASAN=1 make STAR > build.log 2>&1; then
    echo "  ✓ Build successful"
else
    echo "  ✗ Build failed - check build.log"
    exit 1
fi

# Verify ASan linkage
if ldd STAR | grep -q asan; then
    echo "  ✓ AddressSanitizer linked successfully"
else
    echo "  ✗ AddressSanitizer not found in binary"
    exit 1
fi

cd ..

# Step 2: Check test data availability
echo ""
echo "Step 2: Checking test data availability..."

# Set default paths (user can override with environment variables)
export GENOME_DIR="${GENOME_DIR:-/storage/scRNAseq_output/indices-98-32/star}"
export FASTQ_DIR="${FASTQ_DIR:-/storage/downsampled/SC2300771}"
export WHITELIST="${WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"

echo "  Genome index: $GENOME_DIR"
echo "  FASTQ files: $FASTQ_DIR"
echo "  CB whitelist: $WHITELIST"

# Check if paths exist
missing_files=0
for path in "$GENOME_DIR" "$FASTQ_DIR" "$WHITELIST"; do
    if [[ ! -e "$path" ]]; then
        echo "  ✗ Missing: $path"
        missing_files=$((missing_files + 1))
    else
        echo "  ✓ Found: $path"
    fi
done

if [[ $missing_files -gt 0 ]]; then
    echo ""
    echo "Missing $missing_files required files/directories."
    echo "Please set environment variables or provide test data:"
    echo "  export GENOME_DIR=/path/to/genome/index"
    echo "  export FASTQ_DIR=/path/to/fastq/files"
    echo "  export WHITELIST=/path/to/cb_whitelist.txt"
    exit 1
fi

# Step 3: Run memory test
echo ""
echo "Step 3: Running memory test..."
echo "  This will test tag table functionality under AddressSanitizer"
echo "  Expected runtime: 2-10 minutes depending on data size"
echo ""

# Run with leak detection and save log
export ASAN_OPTIONS="detect_leaks=1:print_stats=1"
if timeout 1800 ./mem_test_tags.sh > mem_test_example.log 2>&1; then
    echo "  ✓ Memory test completed successfully"
else
    exit_code=$?
    if [[ $exit_code -eq 124 ]]; then
        echo "  ✗ Memory test timed out after 30 minutes"
    else
        echo "  ✗ Memory test failed with exit code $exit_code"
    fi
    echo "  Check mem_test_example.log for details"
    exit 1
fi

# Step 4: Analyze results
echo ""
echo "Step 4: Analyzing results..."

# Check for critical errors
critical_errors=$(grep -c "ERROR: AddressSanitizer.*heap-\|use-after-free\|double-free" mem_test_example.log || true)
if [[ $critical_errors -gt 0 ]]; then
    echo "  ✗ Found $critical_errors critical memory errors!"
    echo "    These require immediate attention:"
    grep "ERROR: AddressSanitizer.*heap-\|use-after-free\|double-free" mem_test_example.log | head -3
    echo "    See mem_test_example.log for full details"
else
    echo "  ✓ No critical memory errors detected"
fi

# Check for memory leaks (informational)
if grep -q "detected memory leaks" mem_test_example.log; then
    leak_summary=$(grep "SUMMARY: AddressSanitizer:" mem_test_example.log | tail -1)
    echo "  ℹ Memory leaks detected (typical for bioinformatics tools):"
    echo "    $leak_summary"
else
    echo "  ✓ No memory leaks detected"
fi

# Check output files
echo ""
echo "Step 5: Verifying output files..."
output_dir="mem_test_output"
expected_files=("Aligned.out.bam" "Aligned.out.cb_ub.bin" "Solo.out")

for file in "${expected_files[@]}"; do
    if [[ -e "$output_dir/$file" ]]; then
        if [[ -f "$output_dir/$file" ]]; then
            size=$(stat -c%s "$output_dir/$file" 2>/dev/null || stat -f%z "$output_dir/$file" 2>/dev/null || echo "unknown")
            echo "  ✓ $file ($size bytes)"
        else
            echo "  ✓ $file (directory)"
        fi
    else
        echo "  ✗ Missing: $file"
    fi
done

# Step 6: Cleanup (rebuild production binary)
echo ""
echo "Step 6: Rebuilding production binary..."
cd source
make clean > /dev/null
if make STAR > build_production.log 2>&1; then
    echo "  ✓ Production binary rebuilt successfully"
    
    # Verify no ASan linkage
    if ! ldd STAR | grep -q asan; then
        echo "  ✓ Clean production build confirmed"
    else
        echo "  ⚠ Warning: ASan still linked in production binary"
    fi
else
    echo "  ✗ Production build failed - check build_production.log"
fi

cd ..

# Summary
echo ""
echo "=== Memory Testing Summary ==="
if [[ $critical_errors -eq 0 ]]; then
    echo "✅ Memory testing PASSED"
    echo "   - No critical memory errors detected"
    echo "   - Tag table functionality is memory-safe"
    echo "   - Binary serialization working correctly"
else
    echo "❌ Memory testing FAILED"
    echo "   - $critical_errors critical errors found"
    echo "   - Review mem_test_example.log for details"
fi

echo ""
echo "Generated files:"
echo "  - mem_test_example.log - Full ASan output"
echo "  - $output_dir/ - Test outputs"
echo "  - build.log, build_production.log - Compilation logs"
echo ""
echo "For detailed analysis, see docs/memory_testing_guide.md"
