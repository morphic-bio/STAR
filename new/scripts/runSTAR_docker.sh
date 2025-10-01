#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: runSTAR_docker.sh \
    --sample-id SAMPLE_ID \
    --base-dir PATH \
    --whitelist PATH \
    --genome-dir PATH \
    --threads N \
    --output-base PATH \
    [--docker-image IMAGE] [--star-cmd CMD] [--local] \
    [--<STAR_FLAG> VALUE ...]

Required arguments:
  --sample-id       Sample identifier (used to locate FASTQs)
  --base-dir        Directory containing FASTQ files
  --whitelist       Cell barcode whitelist file
  --genome-dir      STAR genome directory
  --threads         Number of threads to allocate to STAR
  --output-base     Parent directory for STAR output

Optional arguments:
  --docker-image    Docker image to run (default: biodepot/star_fork:latest)
  --star-cmd        STAR executable name (default: STAR_fork)
  --local           Run locally instead of via Docker
  -h, --help        Show this help message and exit

Additional STAR flags:
  Any unrecognized --flag will be passed to STAR. If the flag already exists
  in the default parameters, its value will be replaced. If it's new, it will
  be appended to the STAR command.
USAGE
}

USE_DOCKER=1
DOCKER_IMAGE="biodepot/star_fork:latest"
STAR_CMD_NAME="STAR"
declare -A EXTRA_FLAGS=()

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-id)
            SAMPLE_ID="$2"
            shift 2
            ;;
        --base-dir)
            BASE_DIR_RAW="$2"
            shift 2
            ;;
        --whitelist)
            WHITELIST_RAW="$2"
            shift 2
            ;;
        --genome-dir)
            GENOME_DIR_RAW="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --output-base)
            OUTPUT_BASE_RAW="$2"
            shift 2
            ;;
        --docker-image)
            DOCKER_IMAGE="$2"
            shift 2
            ;;
        --star-cmd)
            STAR_CMD_NAME="$2"
            shift 2
            ;;
        --local)
            USE_DOCKER=0
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --*)
            # Handle extra STAR flags
            flag="$1"
            shift
            values=()
            # Collect all values until next -- flag or end of arguments
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                values+=("$1")
                shift
            done
            EXTRA_FLAGS["$flag"]="${values[*]}"
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

for var_name in SAMPLE_ID BASE_DIR_RAW WHITELIST_RAW GENOME_DIR_RAW THREADS OUTPUT_BASE_RAW; do
    if [[ -z "${!var_name:-}" ]]; then
        echo "Missing required argument: ${var_name/_RAW/}" >&2
        usage >&2
        exit 1
    fi
done

if ! [[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]]; then
    echo "--threads must be a positive integer" >&2
    exit 1
fi

resolve_path() {
    local target="$1"
    if command -v realpath >/dev/null 2>&1; then
        realpath "$target"
    elif command -v python3 >/dev/null 2>&1; then
        python3 -c 'import os, sys; print(os.path.realpath(sys.argv[1]))' "$target"
    else
        local dir
        dir=$(cd "$(dirname "$target")" && pwd)
        echo "$dir/$(basename "$target")"
    fi
}

if [[ ! -d "$BASE_DIR_RAW" ]]; then
    echo "Input directory not found: $BASE_DIR_RAW" >&2
    exit 1
fi
BASE_DIR=$(resolve_path "$BASE_DIR_RAW")

if [[ ! -f "$WHITELIST_RAW" ]]; then
    echo "Whitelist not found: $WHITELIST_RAW" >&2
    exit 1
fi
WHITELIST=$(resolve_path "$WHITELIST_RAW")

if [[ ! -d "$GENOME_DIR_RAW" ]]; then
    echo "Genome directory not found: $GENOME_DIR_RAW" >&2
    exit 1
fi
GENOME_DIR=$(resolve_path "$GENOME_DIR_RAW")

if [[ ! -d "$OUTPUT_BASE_RAW" ]]; then
    echo "Output base directory not found: $OUTPUT_BASE_RAW" >&2
    exit 1
fi
OUTPUT_BASE=$(resolve_path "$OUTPUT_BASE_RAW")

NEW_DIR="${OUTPUT_BASE}/${SAMPLE_ID}"
TEMP_DIR="${NEW_DIR}/tmp"

ulimit -n 100000 || true

cleanup_dir() {
    local dir="$1"
    echo "Cleaning up: $dir"
    rm -rf "$dir" "${dir}_STARtmp" 2>/dev/null || true
    mkdir -p "$dir"
}

echo "=== Checking inputs/binaries ==="
if (( USE_DOCKER )); then
    if ! command -v docker >/dev/null 2>&1; then
        echo "Docker is not available but required unless --local is specified" >&2
        exit 1
    fi
else
    if ! command -v "$STAR_CMD_NAME" >/dev/null 2>&1; then
        echo "STAR executable not found in PATH: $STAR_CMD_NAME" >&2
        exit 1
    fi
fi

declare -a R1_LIST=()
declare -a R2_LIST=()
shopt -s nullglob
# Find all R2 files matching the pattern
while IFS= read -r -d '' file; do
    R2_LIST+=("$file")
done < <(find "${BASE_DIR}" -name "${SAMPLE_ID}_*_R2_001.fastq.gz" -print0)

# Create R1_LIST by converting R2 filenames to R1 and checking existence
for r2_file in "${R2_LIST[@]}"; do
    r1_file="${r2_file/_R2_/_R1_}"
    if [[ -f "$r1_file" ]]; then
        R1_LIST+=("$r1_file")
    else
        echo "WARNING: Corresponding R1 file not found for $r2_file" >&2
        #exit if the R1 file is not found
        exit 1
    fi
done
shopt -u nullglob

if (( ${#R2_LIST[@]} == 0 )) || (( ${#R1_LIST[@]} == 0 )); then
    echo "ERROR: Could not find R1/R2 FASTQs in $BASE_DIR for sample $SAMPLE_ID" >&2
    exit 1
fi

join_by() {
    local IFS="$1"
    shift
    echo "$*"
}

R2_FILES=$(join_by ',' "${R2_LIST[@]}")
R1_FILES=$(join_by ',' "${R1_LIST[@]}")

echo "R2: $R2_FILES"
echo "R1: $R1_FILES"

COMMON_PARAMS=(
    --runThreadN "$THREADS"
    --outTmpDir "$TEMP_DIR"
    --soloType CB_UMI_Simple
    --soloCBlen 16
    --soloUMIlen 12
    --soloUMIstart 17
    --soloCBstart 1
    --soloBarcodeReadLength 0
    --soloCBwhitelist "$WHITELIST"
    --genomeDir "$GENOME_DIR"
    --limitIObufferSize 50000000 50000000
    --outSJtype None
    --outBAMcompression 6
    --soloMultiMappers Unique
    --alignIntronMax 1
    --alignMatesGapMax 0
    --outFilterMismatchNmax 10
    --outFilterMismatchNoverReadLmax 1.0
    --outFilterMatchNmin 16
    --outSAMunmapped None
    --outFilterMatchNminOverLread 0
    --outFilterMultimapNmax 10000
    --outFilterMultimapScoreRange 1000
    --outSAMmultNmax 10000
    --winAnchorMultimapNmax 200
    --outSAMprimaryFlag AllBestScore
    --outFilterScoreMin 0
    --outFilterScoreMinOverLread 0
    --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN
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
)

# Process extra flags
if [[ ${#EXTRA_FLAGS[@]} -gt 0 ]]; then
    echo "=== Processing extra STAR flags ==="
    
    # Create a temporary array to build the final parameters
    declare -a FINAL_PARAMS=()
    declare -A PROCESSED_FLAGS=()
    
    # First pass: copy existing parameters and mark which flags we've seen
    i=0
    while [[ $i -lt ${#NEW_TAG_TABLE_PARAMS[@]} ]]; do
        param="${NEW_TAG_TABLE_PARAMS[$i]}"
        if [[ "$param" =~ ^-- ]]; then
            # This is a flag
            if [[ -n "${EXTRA_FLAGS[$param]:-}" ]]; then
                # We have a replacement for this flag
                extra_values="${EXTRA_FLAGS[$param]}"
                if [[ $((i+1)) -lt ${#NEW_TAG_TABLE_PARAMS[@]} && ! "${NEW_TAG_TABLE_PARAMS[$((i+1))]}" =~ ^-- ]]; then
                    # Original has a value, check if it's the same
                    original_value="${NEW_TAG_TABLE_PARAMS[$((i+1))]}"
                    if [[ "$extra_values" == "$original_value" ]]; then
                        echo "Flag $param: identical value '$extra_values', keeping original"
                        FINAL_PARAMS+=("$param" "$original_value")
                        i=$((i+2))
                    else
                        echo "Flag $param: replacing '$original_value' with '$extra_values'"
                        FINAL_PARAMS+=("$param")
                        if [[ -n "$extra_values" ]]; then
                            # Split extra_values by spaces and add each as separate parameter
                            read -ra values_array <<< "$extra_values"
                            FINAL_PARAMS+=("${values_array[@]}")
                        fi
                        i=$((i+2))
                    fi
                else
                    # Original has no value
                    echo "Flag $param: adding value '$extra_values'"
                    FINAL_PARAMS+=("$param")
                    if [[ -n "$extra_values" ]]; then
                        read -ra values_array <<< "$extra_values"
                        FINAL_PARAMS+=("${values_array[@]}")
                    fi
                    i=$((i+1))
                fi
                PROCESSED_FLAGS["$param"]=1
            else
                # Keep original parameter
                FINAL_PARAMS+=("$param")
                if [[ $((i+1)) -lt ${#NEW_TAG_TABLE_PARAMS[@]} && ! "${NEW_TAG_TABLE_PARAMS[$((i+1))]}" =~ ^-- ]]; then
                    FINAL_PARAMS+=("${NEW_TAG_TABLE_PARAMS[$((i+1))]}")
                    i=$((i+2))
                else
                    i=$((i+1))
                fi
            fi
        else
            # This shouldn't happen in our parameter structure, but handle it
            FINAL_PARAMS+=("$param")
            i=$((i+1))
        fi
    done
    
    # Second pass: add any extra flags that weren't already present
    for flag in "${!EXTRA_FLAGS[@]}"; do
        if [[ -z "${PROCESSED_FLAGS[$flag]:-}" ]]; then
            echo "Flag $flag: adding new flag with value '${EXTRA_FLAGS[$flag]}'"
            FINAL_PARAMS+=("$flag")
            if [[ -n "${EXTRA_FLAGS[$flag]}" ]]; then
                read -ra values_array <<< "${EXTRA_FLAGS[$flag]}"
                FINAL_PARAMS+=("${values_array[@]}")
            fi
        fi
    done
    
    NEW_TAG_TABLE_PARAMS=("${FINAL_PARAMS[@]}")
fi

cleanup_dir "$NEW_DIR"
mkdir -p "$TEMP_DIR"

if (( USE_DOCKER )); then
    echo "=== Running STAR via Docker image ${DOCKER_IMAGE} ==="
    declare -A MOUNT_MODES=()
    add_mount() {
        local path="$1"
        local mode="$2"
        [[ -z "$path" ]] && return
        if [[ -z "${MOUNT_MODES[$path]:-}" || "$mode" == "rw" ]]; then
            MOUNT_MODES[$path]="$mode"
        fi
    }

    add_mount "$BASE_DIR" ro
    add_mount "$(dirname "$WHITELIST")" ro
    add_mount "$GENOME_DIR" ro
    add_mount "$OUTPUT_BASE" rw

    DOCKER_CMD=(docker run --rm)
    DOCKER_CMD+=(-u "$(id -u)":"$(id -g)")
    for path in "${!MOUNT_MODES[@]}"; do
        DOCKER_CMD+=(-v "$path:$path:${MOUNT_MODES[$path]}")
    done
    DOCKER_CMD+=("$DOCKER_IMAGE" /bin/bash -c)
    DOCKER_CMD+=("rm -rf $TEMP_DIR && $STAR_CMD_NAME \"\$@\"" --)
    echo "${DOCKER_CMD[@]}"
    "${DOCKER_CMD[@]}" \
        "${NEW_TAG_TABLE_PARAMS[@]}" \
        --outFileNamePrefix "${NEW_DIR}/"
else
    echo "=== Running STAR locally with ${STAR_CMD_NAME} ==="
    "$STAR_CMD_NAME" \
        "${NEW_TAG_TABLE_PARAMS[@]}" \
        --outFileNamePrefix "${NEW_DIR}/"
fi

echo "Run finished"
