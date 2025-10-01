#!/usr/bin/env bash
set -euo pipefail

# Minimal test script for tag-table functionality with AddressSanitizer
# Based on mem_plan.txt Phase 3
#
# For detailed usage instructions, see:
#   docs/memory_testing_guide.md - Complete documentation
#   MEMORY_TESTING.md - Quick reference
#
# Prerequisites:
#   1. Build STAR with ASan: cd source && ASAN=1 make clean && ASAN=1 make STAR
#   2. Set environment variables (see documentation for defaults)
#   3. Run: ASAN_OPTIONS="detect_leaks=1" ./mem_test_tags.sh

# Configuration with defaults
STAR_BIN=${STAR_BIN:-./source/STAR}
GENOME_DIR=${GENOME_DIR:-/storage/scRNAseq_output/indices-98-32/star}
FASTQ_DIR=${FASTQ_DIR:-/storage/downsampled/SC2300771}
WHITELIST=${WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}
OUTDIR=${OUTDIR:-mem_test_output}
THREADS=${THREADS:-4}

echo "=== Memory Test for Tag Table Functionality ==="
echo "STAR Binary: $STAR_BIN"
echo "Genome Dir: $GENOME_DIR"
echo "FASTQ Dir: $FASTQ_DIR"
echo "Whitelist: $WHITELIST"
echo "Output Dir: $OUTDIR"
echo "Threads: $THREADS"

# Check requirements
[[ -f "$STAR_BIN" ]] || { echo "ERROR: STAR binary not found: $STAR_BIN"; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: Genome dir not found: $GENOME_DIR"; exit 1; }
[[ -d "$FASTQ_DIR" ]] || { echo "ERROR: FASTQ dir not found: $FASTQ_DIR"; exit 1; }
[[ -f "$WHITELIST" ]] || { echo "ERROR: Whitelist not found: $WHITELIST"; exit 1; }

# Build FASTQ file lists (using first 2 lanes for speed)
R2_FILES=""
R1_FILES=""
for lane in L001 L002; do
    r2_pattern="${FASTQ_DIR}/SC2300771_*_${lane}_R2_001.fastq.gz"
    r1_pattern="${FASTQ_DIR}/SC2300771_*_${lane}_R1_001.fastq.gz"
    
    if ls $r2_pattern 1>/dev/null 2>&1; then
        [[ -n "$R2_FILES" ]] && R2_FILES+=","
        R2_FILES+="$(ls $r2_pattern)"
    fi
    
    if ls $r1_pattern 1>/dev/null 2>&1; then
        [[ -n "$R1_FILES" ]] && R1_FILES+=","
        R1_FILES+="$(ls $r1_pattern)"
    fi
done

[[ -n "$R2_FILES" && -n "$R1_FILES" ]] || { echo "ERROR: Could not find R1/R2 FASTQs in $FASTQ_DIR"; exit 1; }

echo "R2 files: $R2_FILES"
echo "R1 files: $R1_FILES"

# Clean up output directory
rm -rf "$OUTDIR" "${OUTDIR}_STARtmp"
mkdir -p "$OUTDIR"

echo ""
echo "=== Running STAR with Tag Table Export (ASan enabled) ==="

# Run STAR with minimal parameters to trigger tag-table code path
"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --genomeDir "$GENOME_DIR" \
  --readFilesIn "$R2_FILES" "$R1_FILES" \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 \
  --soloUMIlen 12 \
  --soloUMIstart 17 \
  --soloCBstart 1 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist "$WHITELIST" \
  --soloMultiMappers Unique \
  --outSAMtype BAM Unsorted \
  --soloAddTagsToUnsorted yes \
  --soloWriteTagTable Default \
  --outSAMattributes NH HI AS nM NM CB UB CR CY UR UY GX GN \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --alignEndsType Local \
  --outFileNamePrefix "$OUTDIR/"

echo ""
echo "=== Memory Test Completed ==="
echo "Output directory: $OUTDIR"
echo "Check for:"
echo "  - Binary tag table: ${OUTDIR}/Aligned.out.cb_ub.bin"
echo "  - BAM with tags: ${OUTDIR}/Aligned.out.bam"
echo "  - Solo matrices: ${OUTDIR}/Solo.out/"

# Quick validation
if [[ -f "${OUTDIR}/Aligned.out.cb_ub.bin" ]]; then
    echo "✓ Binary tag table created successfully"
    ls -lh "${OUTDIR}/Aligned.out.cb_ub.bin"
else
    echo "✗ Binary tag table not found"
    exit 1
fi

if [[ -f "${OUTDIR}/Aligned.out.bam" ]]; then
    echo "✓ BAM file created successfully"
    ls -lh "${OUTDIR}/Aligned.out.bam"
else
    echo "✗ BAM file not found"
    exit 1
fi

echo ""
echo "Memory test completed successfully!"
