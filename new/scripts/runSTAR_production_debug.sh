#!/usr/bin/env bash

set -euo pipefail

# User-editable parameters
STAR_BIN="/mnt/pikachu/STAR/source/STAR"
THREADS=24
SAMPLE_ID="SC2300771"
WHITELIST="/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt"
GENOME_DIR="/storage/scRNAseq_output/indices-98-32/star"
BASE_DIR="/storage/JAX_sequences/${SAMPLE_ID}"
OUT_PREFIX="/storage/Alignments/${SAMPLE_ID}_debug/"
TEMP_DIR="/storage/tmp/${SAMPLE_ID}_debug"

# Find all R1 files and verify corresponding R2 files exist
R1_FILES=""
R2_FILES=""

while IFS= read -r -d '' r1_file; do
    # Convert R1 filename to R2 filename
    r2_file="${r1_file/_R1_/_R2_}"
    
    # Check if corresponding R2 file exists
    if [[ ! -f "$r2_file" ]]; then
        echo "ERROR: Found R1 file but missing corresponding R2 file:"
        echo "  R1: $r1_file"
        echo "  R2: $r2_file (missing)"
        exit 1
    fi
    
    # Add to file lists
    [[ -n "$R1_FILES" ]] && R1_FILES+=","
    R1_FILES+="$r1_file"
    [[ -n "$R2_FILES" ]] && R2_FILES+=","
    R2_FILES+="$r2_file"
    
done < <(find "$BASE_DIR" -name "${SAMPLE_ID}_*_R1_001.fastq.gz" -print0 | sort -z)

[[ -n "$R1_FILES" && -n "$R2_FILES" ]] || { echo "ERROR: Could not find any R1/R2 FASTQs in $BASE_DIR"; exit 1; }

# Count files found
R1_COUNT=$(echo "$R1_FILES" | tr ',' '\n' | wc -l)
R2_COUNT=$(echo "$R2_FILES" | tr ',' '\n' | wc -l)

echo "Found $R1_COUNT R1 files and $R2_COUNT R2 files"
echo "R1: $R1_FILES"
echo "R2: $R2_FILES"
rm -rf "$TEMP_DIR"

ulimit -n 100000 || true


"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --outTmpDir "$TEMP_DIR" \
  --readMapNumber -1 \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 \
  --soloUMIlen 12 \
  --soloUMIstart 17 \
  --soloCBstart 1 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist "$WHITELIST" \
  --genomeDir "$GENOME_DIR" \
  --limitIObufferSize 25000000 25000000 \
  --outSJtype None \
  --outBAMcompression 6 \
  --soloMultiMappers Unique \
  --alignIntronMax 1 \
  --alignMatesGapMax 0 \
  --outFilterMismatchNmax 6 \
  --outFilterMismatchNoverReadLmax 1.0 \
  --outFilterMatchNmin 25 \
  --outSAMunmapped None \
  --outFilterMatchNminOverLread 0 \
  --outFilterMultimapNmax 10000 \
  --outFilterMultimapScoreRange 1 \
  --outSAMmultNmax 10000 \
  --winAnchorMultimapNmax 200 \
  --outSAMprimaryFlag AllBestScore \
  --outFilterScoreMin 0 \
  --outFilterScoreMinOverLread 0 \
  --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --readFilesIn "$R2_FILES" "$R1_FILES" \
  --alignEndsType Local \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --soloWriteTagTable Default \
  --soloAddTagsToUnsorted no \
  --outFileNamePrefix "$OUT_PREFIX"


