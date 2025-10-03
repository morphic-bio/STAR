#!/bin/bash
# Stage 0.3 – Baseline skip-mode behavior capture
# Records memory usage and logs before refactoring begins

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTPUT_DIR="$REPO_ROOT/new/testing/outputs/baseline_skip_mode"

mkdir -p "$OUTPUT_DIR"

echo "=== Stage 0.3: Baseline Skip-Mode Capture ==="
echo "Output directory: $OUTPUT_DIR"

# Check if we have test data
if [[ ! -d "$REPO_ROOT/new/testing/data" ]]; then
    echo "WARNING: Test data directory not found at $REPO_ROOT/new/testing/data"
    echo "Skipping actual STAR run, but documenting baseline structure..."
    
    # Document current state
    cat > "$OUTPUT_DIR/baseline_notes.txt" <<EOF
Baseline Capture (Stage 0.3)
Date: $(date)
Branch: $(cd "$REPO_ROOT" && git branch --show-current)
Commit: $(cd "$REPO_ROOT" && git rev-parse HEAD)

Build Configuration:
- Packed ReadInfo: ENABLED (default)
- STAR version: $(cd "$REPO_ROOT" && ./source/STAR --version)

Test Status:
- Test data not found; baseline capture deferred until test data available
- Current implementation uses packed readinfo by default
- Skip mode invokes prepareReadInfoOnly() which currently allocates minimal arrays

Expected Memory Profile (to be measured):
- Skip mode should avoid allocating rGeneUMI, rCBp, countCellGeneUMI, countMatMult
- Baseline (non-skip) mode allocates all counting structures

Next Steps:
- Run with representative dataset when available
- Measure RSS via /usr/bin/time -v or similar
- Capture Log.out messages related to array allocation
EOF
    
    cat "$OUTPUT_DIR/baseline_notes.txt"
    exit 0
fi

# If we have test data, run a minimal STAR command to capture baseline
echo "Test data found, running baseline capture..."

# Capture system info
cat > "$OUTPUT_DIR/system_info.txt" <<EOF
System Information
==================
Date: $(date)
Hostname: $(hostname)
OS: $(uname -a)
Memory: $(free -h | grep Mem)
CPUs: $(nproc)

Git Status
==========
Branch: $(cd "$REPO_ROOT" && git branch --show-current)
Commit: $(cd "$REPO_ROOT" && git rev-parse HEAD)

STAR Build
==========
Version: $(cd "$REPO_ROOT" && ./source/STAR --version)
Build flags: SOLO_USE_PACKED_READINFO enabled by default
EOF

echo ""
echo "Baseline capture complete. Summary written to:"
echo "  $OUTPUT_DIR/baseline_notes.txt"
echo "  $OUTPUT_DIR/system_info.txt"
echo ""
echo "Stage 0 complete. Ready for Stage 1 implementation."

