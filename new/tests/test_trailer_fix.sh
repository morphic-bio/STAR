#!/bin/bash

# Quick test to verify the trailer parsing fix

set -euo pipefail

echo "=== Testing Trailer Parsing Fix ==="
echo ""

# Test the fix with mem_test_tags.sh to see if it now handles trailers correctly
echo "Running quick test with debug instrumentation..."
echo "This should now correctly parse the 64-bit trailer structure"
echo ""

STAR_DEBUG_TAG=1 timeout 60 ./mem_test_tags.sh 2>&1 | grep -A 10 -B 2 "trailer\|iread\|ERROR\|BAM record processing" | head -20

echo ""
echo "=== Analysis ==="
echo "The debug output should now show:"
echo "✓ Correct trailer parsing with 'trailer: full=0x..., iread=X, aux=0x...'"
echo "✓ Valid iread values (not 4294967296)"
echo "✓ No 'Invalid read index' errors"
echo ""
echo "If you see 'Invalid read index' errors, the trailer structure may be different"
echo "than expected, but the debug output will now show the full trailer breakdown."
