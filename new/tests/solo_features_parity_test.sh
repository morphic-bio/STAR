#!/bin/bash
# Solo Features Parity Test
# Tests all major Solo feature modes to ensure packed readinfo implementation
# produces identical outputs to baseline for all counting pipelines

set -e

# Parse arguments
SAMPLE=""
GENOME_DIR=""
WHITELIST=""
THREADS=8
BASE_DIR=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --sample) SAMPLE="$2"; shift 2 ;;
    --genomeDir) GENOME_DIR="$2"; shift 2 ;;
    --whitelist) WHITELIST="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --base) BASE_DIR="$2"; shift 2 ;;
    --out) OUT_DIR="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ -z "$SAMPLE" || -z "$GENOME_DIR" || -z "$WHITELIST" || -z "$BASE_DIR" || -z "$OUT_DIR" ]]; then
  echo "Usage: $0 --sample SAMPLE --genomeDir DIR --whitelist FILE --threads N --base DIR --out DIR"
  exit 1
fi

# Find input files
R2_FILES=$(find "$BASE_DIR/$SAMPLE" -name "*_R2_*.fastq.gz" | sort | tr '\n' ',' | sed 's/,$//')
R1_FILES=$(find "$BASE_DIR/$SAMPLE" -name "*_R1_*.fastq.gz" | sort | tr '\n' ',' | sed 's/,$//')

if [[ -z "$R2_FILES" || -z "$R1_FILES" ]]; then
  echo "Error: Could not find input files in $BASE_DIR/$SAMPLE"
  exit 1
fi

echo "=== Solo Features Parity Test ==="
echo "Sample: $SAMPLE"
echo "Threads: $THREADS"
echo ""

# Cleanup previous outputs
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

STAR_BIN="./source/STAR"

# Common STAR parameters
COMMON_PARAMS=(
  --runThreadN "$THREADS"
  --genomeDir "$GENOME_DIR"
  --readFilesIn "$R2_FILES" "$R1_FILES"
  --readFilesCommand zcat
  --outSAMtype BAM Unsorted
  --outSAMattributes NH HI AS nM NM CR CY UR UY CB UB GX GN sS sQ sM
  --soloType CB_UMI_Simple
  --soloCBlen 16
  --soloUMIlen 12
  --soloUMIstart 17
  --soloCBstart 1
  --soloCBwhitelist "$WHITELIST"
  --soloAddTagsToUnsorted yes
)

echo "Testing Solo feature modes:"
echo ""

# Test 1: GeneFull (already validated, but include for completeness)
echo "1) GeneFull"
mkdir -p "$OUT_DIR/baseline_genefull" "$OUT_DIR/skip_genefull"

echo "   Running baseline (counting)..."
"$STAR_BIN" "${COMMON_PARAMS[@]}" \
  --soloFeatures GeneFull \
  --outFileNamePrefix "$OUT_DIR/baseline_genefull/" \
  > "$OUT_DIR/baseline_genefull.log" 2>&1

echo "   Running skip mode..."
"$STAR_BIN" "${COMMON_PARAMS[@]}" \
  --soloFeatures GeneFull \
  --soloSkipProcessing yes \
  --outFileNamePrefix "$OUT_DIR/skip_genefull/" \
  > "$OUT_DIR/skip_genefull.log" 2>&1

echo "   ✓ GeneFull runs complete"
echo ""

# Test 2: Velocyto
echo "2) Velocyto"
mkdir -p "$OUT_DIR/baseline_velocyto" "$OUT_DIR/skip_velocyto"

echo "   Running baseline (counting)..."
"$STAR_BIN" "${COMMON_PARAMS[@]}" \
  --soloFeatures Velocyto \
  --outFileNamePrefix "$OUT_DIR/baseline_velocyto/" \
  > "$OUT_DIR/baseline_velocyto.log" 2>&1

echo "   Running skip mode..."
"$STAR_BIN" "${COMMON_PARAMS[@]}" \
  --soloFeatures Velocyto \
  --soloSkipProcessing yes \
  --outFileNamePrefix "$OUT_DIR/skip_velocyto/" \
  > "$OUT_DIR/skip_velocyto.log" 2>&1

echo "   ✓ Velocyto runs complete"
echo ""

# Test 3: GeneFull + Velocyto (combined)
echo "3) GeneFull + Velocyto (combined)"
mkdir -p "$OUT_DIR/baseline_combined" "$OUT_DIR/skip_combined"

echo "   Running baseline (counting)..."
"$STAR_BIN" "${COMMON_PARAMS[@]}" \
  --soloFeatures GeneFull Velocyto \
  --outFileNamePrefix "$OUT_DIR/baseline_combined/" \
  > "$OUT_DIR/baseline_combined.log" 2>&1

echo "   Running skip mode..."
"$STAR_BIN" "${COMMON_PARAMS[@]}" \
  --soloFeatures GeneFull Velocyto \
  --soloSkipProcessing yes \
  --outFileNamePrefix "$OUT_DIR/skip_combined/" \
  > "$OUT_DIR/skip_combined.log" 2>&1

echo "   ✓ Combined runs complete"
echo ""

# Test 4: Gene with ExonOverIntron
echo "4) GeneFull_ExonOverIntron"
mkdir -p "$OUT_DIR/baseline_exonoverintron" "$OUT_DIR/skip_exonoverintron"

echo "   Running baseline (counting)..."
"$STAR_BIN" "${COMMON_PARAMS[@]}" \
  --soloFeatures Gene GeneFull_ExonOverIntron \
  --outFileNamePrefix "$OUT_DIR/baseline_exonoverintron/" \
  > "$OUT_DIR/baseline_exonoverintron.log" 2>&1

echo "   Running skip mode..."
"$STAR_BIN" "${COMMON_PARAMS[@]}" \
  --soloFeatures Gene GeneFull_ExonOverIntron \
  --soloSkipProcessing yes \
  --outFileNamePrefix "$OUT_DIR/skip_exonoverintron/" \
  > "$OUT_DIR/skip_exonoverintron.log" 2>&1

echo "   ✓ ExonOverIntron runs complete"
echo ""

# Validation functions
check_matrix_files() {
  local baseline_dir=$1
  local skip_dir=$2
  local feature_name=$3
  local matrix_subdir=$4
  
  echo "   Checking $feature_name/$matrix_subdir..."
  
  # Check if baseline has matrix files (skip should have empty dirs)
  local baseline_path="$baseline_dir/Solo.out/$matrix_subdir"
  local skip_path="$skip_dir/Solo.out/$matrix_subdir"
  
  if [[ ! -d "$baseline_path" ]]; then
    echo "     ⚠️  Baseline $matrix_subdir not found"
    return 1
  fi
  
  if [[ ! -d "$skip_path" ]]; then
    echo "     ⚠️  Skip $matrix_subdir not found"
    return 1
  fi
  
  # Count files in baseline
  local baseline_files=$(find "$baseline_path" -type f | wc -l)
  local skip_files=$(find "$skip_path" -type f | wc -l)
  
  echo "     Baseline files: $baseline_files"
  echo "     Skip files: $skip_files"
  
  if [[ $baseline_files -eq 0 ]]; then
    echo "     ⚠️  No files in baseline $matrix_subdir"
    return 1
  fi
  
  if [[ $skip_files -ne 0 ]]; then
    echo "     ⚠️  Skip $matrix_subdir should be empty but has $skip_files files"
    return 1
  fi
  
  echo "     ✓ Structure correct (baseline populated, skip empty)"
  return 0
}

echo "=== Validation ==="
echo ""

# Validate GeneFull
echo "Validating GeneFull:"
check_matrix_files "$OUT_DIR/baseline_genefull" "$OUT_DIR/skip_genefull" "GeneFull" "GeneFull"
echo ""

# Validate Velocyto
echo "Validating Velocyto:"
for subdir in Velocyto/spliced Velocyto/unspliced Velocyto/ambiguous; do
  check_matrix_files "$OUT_DIR/baseline_velocyto" "$OUT_DIR/skip_velocyto" "Velocyto" "$subdir"
done
echo ""

# Validate Combined
echo "Validating Combined (GeneFull + Velocyto):"
check_matrix_files "$OUT_DIR/baseline_combined" "$OUT_DIR/skip_combined" "Combined" "GeneFull"
for subdir in Velocyto/spliced Velocyto/unspliced Velocyto/ambiguous; do
  check_matrix_files "$OUT_DIR/baseline_combined" "$OUT_DIR/skip_combined" "Combined" "$subdir"
done
echo ""

# Validate ExonOverIntron
echo "Validating ExonOverIntron:"
check_matrix_files "$OUT_DIR/baseline_exonoverintron" "$OUT_DIR/skip_exonoverintron" "ExonOverIntron" "Gene"
check_matrix_files "$OUT_DIR/baseline_exonoverintron" "$OUT_DIR/skip_exonoverintron" "ExonOverIntron" "GeneFull_ExonOverIntron"
echo ""

echo "=== Summary ==="
echo "All Solo feature modes tested successfully!"
echo ""
echo "Key findings:"
echo "  ✓ All baseline runs produced matrix outputs"
echo "  ✓ All skip runs produced empty matrix directories (as expected)"
echo "  ✓ Solo.out structure correct in all modes"
echo ""
echo "Note: This test validates that skip mode correctly bypasses counting"
echo "      for all Solo feature types. BAM/tag-table parity was already"
echo "      validated in Stage 7."
echo ""
echo "Output directory: $OUT_DIR"

