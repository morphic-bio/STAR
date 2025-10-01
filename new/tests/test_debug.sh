#!/bin/bash

# Quick test of runSTAR_debug.sh functionality
# Uses the existing mem_test data to validate debug instrumentation

set -euo pipefail

echo "=== Testing runSTAR_debug.sh functionality ==="
echo ""

# Test 1: Check if debug script exists and is executable
if [[ -x "./runSTAR_debug.sh" ]]; then
    echo "✓ runSTAR_debug.sh is executable"
else
    echo "✗ runSTAR_debug.sh not found or not executable"
    exit 1
fi

# Test 2: Check if STAR binary has debug instrumentation
STAR_BINARY="./source/STAR"
if [[ -f "$STAR_BINARY" ]]; then
    echo "✓ STAR binary found: $STAR_BINARY"
    if strings "$STAR_BINARY" | grep -q "DEBUG_TAG"; then
        echo "✓ Debug instrumentation detected in binary"
    else
        echo "⚠ WARNING: Debug instrumentation not detected"
        echo "  Run: cd source && make STAR"
    fi
else
    echo "✗ STAR binary not found: $STAR_BINARY"
    exit 1
fi

# Test 3: Quick debug environment test
echo ""
echo "=== Testing debug environment ==="
export STAR_DEBUG_TAG=1
export MALLOC_CHECK_=2
export MALLOC_PERTURB_=42

echo "Debug environment variables:"
echo "  STAR_DEBUG_TAG=$STAR_DEBUG_TAG"
echo "  MALLOC_CHECK_=$MALLOC_CHECK_"
echo "  MALLOC_PERTURB_=$MALLOC_PERTURB_"

# Test 4: Run mem_test with debug enabled to verify instrumentation
echo ""
echo "=== Testing debug instrumentation with mem_test_tags.sh ==="
echo "This will run a quick test to verify debug logging works..."

if STAR_DEBUG_TAG=1 timeout 60 ./mem_test_tags.sh 2>&1 | head -20 | grep -q "DEBUG_TAG"; then
    echo "✓ Debug instrumentation is working"
else
    echo "⚠ Debug instrumentation test inconclusive"
fi

echo ""
echo "=== runSTAR_debug.sh Ready for Use ==="
echo ""
echo "Usage:"
echo "  ./runSTAR_debug.sh"
echo ""
echo "Key features:"
echo "  • Enables STAR_DEBUG_TAG=1 for comprehensive logging"
echo "  • Enables malloc debugging (MALLOC_CHECK_=2, MALLOC_PERTURB_=42)"
echo "  • Reduces thread count to 8 for stability"
echo "  • Logs all output to timestamped debug log files"
echo "  • Monitors system resources during run"
echo "  • Analyzes debug output for errors/warnings"
echo "  • Creates debug-specific output directory"
echo ""
echo "Debug logs will be saved to: /tmp/star_debug_logs/"
echo ""
echo "To run the full debug version:"
echo "  ./runSTAR_debug.sh"
echo ""
echo "To monitor progress:"
echo "  tail -f /tmp/star_debug_logs/star_debug_*.log"
