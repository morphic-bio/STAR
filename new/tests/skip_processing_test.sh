#!/usr/bin/env bash

set -euo pipefail

# Skip-processing validation:
# - Runs a baseline STARsolo (full processing)
# - Runs with --soloSkipProcessing yes (minimal processing)
# - Confirms per-read CB/UB outputs match and Solo matrices are omitted in skip mode

usage() {
    cat <<EOF
Usage: $(basename "$0") \
  --sample SC2300771 \
  --genomeDir /storage/flex_filtered_reference/star_index \
  --whitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
  [--threads 24] [--base /storage/downsampled] [--out /storage/downsampled/results]

Notes:
 - Gathers ALL R1/R2 FASTQs at: <base>/<sample> (supports *R1*/*R2* and *_1*/*_2* patterns, .fastq[.gz]|.fq[.gz])
 - Uses the STAR binary from this repo: source/STAR
 - Produces two runs: baseline and skip
EOF
}

SAMPLE=""
GENOME_DIR=""
WHITELIST=""
THREADS=24
BASE=/storage/downsampled
OUT_BASE=/storage/downsampled/results

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --genomeDir) GENOME_DIR="$2"; shift 2;;
    --whitelist) WHITELIST="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --base) BASE="$2"; shift 2;;
    --out) OUT_BASE="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown argument: $1"; usage; exit 1;;
  esac
done

[[ -n "$SAMPLE" && -n "$GENOME_DIR" && -n "$WHITELIST" ]] || { usage; exit 1; }

# Resolve project root and binaries
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../.. && pwd)"
STAR_BIN="${ROOT_DIR}/source/STAR"
DECODER="${ROOT_DIR}/tools/decode_tag_binary"

[[ -x "$STAR_BIN" ]] || { echo "ERROR: STAR binary not found at $STAR_BIN"; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found in PATH"; exit 1; }

# Inputs
IN_DIR="${BASE}/${SAMPLE}"
[[ -d "$IN_DIR" ]] || { echo "ERROR: sample directory not found: $IN_DIR"; exit 1; }

# Gather all R1/R2 files (R2 = cDNA, R1 = barcode for 10x)
mapfile -t R1_ARR < <(find "$IN_DIR" -maxdepth 1 -type f \
  \( -iname '*R1*.fastq.gz' -o -iname '*R1*.fq.gz' -o -iname '*_1.fastq.gz' -o -iname '*_1.fq.gz' \
     -o -iname '*R1*.fastq'    -o -iname '*R1*.fq'    -o -iname '*_1.fastq'    -o -iname '*_1.fq' \) | sort)
mapfile -t R2_ARR < <(find "$IN_DIR" -maxdepth 1 -type f \
  \( -iname '*R2*.fastq.gz' -o -iname '*R2*.fq.gz' -o -iname '*_2.fastq.gz' -o -iname '*_2.fq.gz' \
     -o -iname '*R2*.fastq'    -o -iname '*R2*.fq'    -o -iname '*_2.fastq'    -o -iname '*_2.fq' \) | sort)

[[ ${#R1_ARR[@]} -gt 0 && ${#R2_ARR[@]} -gt 0 ]] || { echo "ERROR: No R1/R2 FASTQs found in $IN_DIR"; exit 1; }
[[ ${#R1_ARR[@]} -eq ${#R2_ARR[@]} ]] || { echo "ERROR: R1/R2 count mismatch: R1=${#R1_ARR[@]} R2=${#R2_ARR[@]}"; printf '%s\n' "R1:" "${R1_ARR[@]}" "R2:" "${R2_ARR[@]}"; exit 1; }

R1_CSV=$(IFS=,; echo "${R1_ARR[*]}")
R2_CSV=$(IFS=,; echo "${R2_ARR[*]}")
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: genomeDir not found: $GENOME_DIR"; exit 1; }
[[ -s "$WHITELIST" ]] || { echo "ERROR: whitelist not found: $WHITELIST"; exit 1; }

# Outputs
BASELINE_OUT="${OUT_BASE}/solo_baseline/"
SKIP_OUT="${OUT_BASE}/solo_skip/"

echo "Cleaning outputs..."
rm -rf "${BASELINE_OUT}" "${BASELINE_OUT}_STARtmp" "${SKIP_OUT}" "${SKIP_OUT}_STARtmp"
mkdir -p "${BASELINE_OUT}" "${SKIP_OUT}"

echo "Building tools/decoder (if needed)..."
make -C "${ROOT_DIR}/tools" >/dev/null
[[ -x "$DECODER" ]] || { echo "ERROR: decoder not found after build: $DECODER"; exit 1; }

run_star() {
  local prefix="$1"
  local skip_flag="$2" # "yes" or "no"
  echo "Running STAR (skip=${skip_flag}) -> ${prefix}"
  "$STAR_BIN" \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$R2_CSV" "$R1_CSV" \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI AS nM NM CR CY UR UY CB UB ZG \
    --soloType CB_UMI_Simple \
    --soloCBlen 16 \
    --soloUMIlen 12 \
    --soloUMIstart 17 \
    --soloCBstart 1 \
    --soloFeatures GeneFull \
    --soloCBwhitelist "$WHITELIST" \
    --soloAddTagsToUnsorted yes \
    --soloWriteTagTable Default \
    --outFileNamePrefix "$prefix" \
    $( [[ "$skip_flag" == "yes" ]] && echo "--soloSkipProcessing yes" )
}

run_star "$BASELINE_OUT" "no"
run_star "$SKIP_OUT" "yes"

echo "Validations:"

echo "1) Solo.out folder structure: baseline vs skip"
[[ -d "${BASELINE_OUT}Solo.out" ]] && echo "  ✓ Baseline Solo.out present" || { echo "  ✗ Baseline Solo.out missing"; exit 1; }
[[ -d "${SKIP_OUT}Solo.out" ]] && echo "  ✓ Skip Solo.out present" || { echo "  ✗ Skip Solo.out missing"; exit 1; }

echo "2) GeneFull folder content check"
# Baseline should have files in GeneFull
BASELINE_GENEFULL="${BASELINE_OUT}Solo.out/GeneFull"
[[ -d "$BASELINE_GENEFULL" ]] || { echo "  ✗ Baseline GeneFull folder missing"; exit 1; }
BASELINE_FILE_COUNT=$(find "$BASELINE_GENEFULL" -type f | wc -l)
if [[ $BASELINE_FILE_COUNT -gt 0 ]]; then
  echo "  ✓ Baseline GeneFull has $BASELINE_FILE_COUNT files"
else
  echo "  ✗ Baseline GeneFull is empty"; exit 1
fi

# Skip should have empty GeneFull folder
SKIP_GENEFULL="${SKIP_OUT}Solo.out/GeneFull"
[[ -d "$SKIP_GENEFULL" ]] || { echo "  ✗ Skip GeneFull folder missing"; exit 1; }
SKIP_FILE_COUNT=$(find "$SKIP_GENEFULL" -type f | wc -l)
if [[ $SKIP_FILE_COUNT -eq 0 ]]; then
  echo "  ✓ Skip GeneFull is empty (no files)"
else
  echo "  ✗ Skip GeneFull has $SKIP_FILE_COUNT files (should be empty)"; exit 1
fi

echo "3) Check BAM files have same number of lines"
BASELINE_BAM_LINES=$(samtools view "${BASELINE_OUT}Aligned.out.bam" | wc -l)
SKIP_BAM_LINES=$(samtools view "${SKIP_OUT}Aligned.out.bam" | wc -l)

if [[ $BASELINE_BAM_LINES -eq $SKIP_BAM_LINES ]]; then
  echo "  ✓ BAM files have same number of lines: $BASELINE_BAM_LINES"
else
  echo "  ✗ BAM line count differs: baseline=$BASELINE_BAM_LINES skip=$SKIP_BAM_LINES"; exit 1
fi

echo "4) Check ZG fields are non-empty in both BAMs"
# Check baseline BAM
BASELINE_ZG_EMPTY=$(samtools view "${BASELINE_OUT}Aligned.out.bam" | \
  awk '{for(i=12;i<=NF;i++){if($i~/^ZG:Z:/){zg=substr($i,6); if(zg=="-" || zg=="") print "empty"}}}' | wc -l)
BASELINE_ZG_NONEMPTY=$(samtools view "${BASELINE_OUT}Aligned.out.bam" | \
  awk '{for(i=12;i<=NF;i++){if($i~/^ZG:Z:/){zg=substr($i,6); if(zg!="-" && zg!="") print "ok"}}}' | wc -l)

# Check skip BAM
SKIP_ZG_EMPTY=$(samtools view "${SKIP_OUT}Aligned.out.bam" | \
  awk '{for(i=12;i<=NF;i++){if($i~/^ZG:Z:/){zg=substr($i,6); if(zg=="-" || zg=="") print "empty"}}}' | wc -l)
SKIP_ZG_NONEMPTY=$(samtools view "${SKIP_OUT}Aligned.out.bam" | \
  awk '{for(i=12;i<=NF;i++){if($i~/^ZG:Z:/){zg=substr($i,6); if(zg!="-" && zg!="") print "ok"}}}' | wc -l)

echo "  Baseline: $BASELINE_ZG_NONEMPTY non-empty ZG, $BASELINE_ZG_EMPTY empty/'-'"
echo "  Skip:     $SKIP_ZG_NONEMPTY non-empty ZG, $SKIP_ZG_EMPTY empty/'-'"

if [[ $BASELINE_ZG_NONEMPTY -gt 0 && $BASELINE_ZG_EMPTY -eq 0 ]]; then
  echo "  ✓ Baseline BAM has non-empty ZG fields"
else
  echo "  ✗ Baseline BAM has empty or missing ZG fields"; exit 1
fi

if [[ $SKIP_ZG_NONEMPTY -gt 0 && $SKIP_ZG_EMPTY -eq 0 ]]; then
  echo "  ✓ Skip BAM has non-empty ZG fields"
else
  echo "  ✗ Skip BAM has empty or missing ZG fields"; exit 1
fi

if [[ $BASELINE_ZG_NONEMPTY -eq $SKIP_ZG_NONEMPTY ]]; then
  echo "  ✓ Both BAMs have same count of non-empty ZG fields"
else
  echo "  ✗ ZG field counts differ between BAMs"; exit 1
fi

echo "5) Compare binary tag-table files with diff (should be identical)"
if diff "${BASELINE_OUT}Aligned.out.cb_ub.bin" "${SKIP_OUT}Aligned.out.cb_ub.bin" >/dev/null; then
  echo "  ✓ Binary tag-table files are identical"
else
  echo "  ✗ Binary tag-table files differ"; exit 1
fi

echo "6) Optional: decode and basic sanity"
"$DECODER" "${BASELINE_OUT}Aligned.out.cb_ub.bin" >/dev/null 2>&1 && echo "  ✓ Baseline tag-table decodes" || echo "  ✗ Baseline tag-table decode failed"
"$DECODER" "${SKIP_OUT}Aligned.out.cb_ub.bin"     >/dev/null 2>&1 && echo "  ✓ Skip tag-table decodes"     || echo "  ✗ Skip tag-table decode failed"

echo "All checks passed!"


