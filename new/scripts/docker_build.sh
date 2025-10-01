#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: docker_build.sh [options]

Build STAR Docker image.

Options:
  --prebuilt              Use pre-built binaries (faster, requires ./build_star.sh first)
  --tag <name>            Docker tag (default: star:latest)
  --no-cache              Build without cache
  --build-args <args>     Additional build arguments
  -h, --help              Show this help

Examples:
  # Option 1: Use pre-built binaries (fastest)
  ./build_star.sh --static --long
  ./docker_build.sh --prebuilt

  # Option 2: Full multi-stage build (slower but self-contained)
  ./docker_build.sh

  # Custom tag
  ./docker_build.sh --tag mystar:v1.0
EOF
}

PREBUILT=0
TAG="biodepot/star_fork:latest"
NO_CACHE=""
BUILD_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prebuilt) PREBUILT=1 ;;
    --tag) TAG="${2:-}"; shift ;;
    --no-cache) NO_CACHE="--no-cache" ;;
    --build-args) BUILD_ARGS="${2:-}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
  shift
done

if (( PREBUILT )); then
  if [[ ! -f "bin/Linux_x86_64_static/STAR" ]]; then
    echo "Error: Pre-built binaries not found. Run './build_star.sh --static --long' first." >&2
    exit 1
  fi
  echo "Building Docker image with pre-built binaries..."
  docker build -f Dockerfile.prebuilt -t "$TAG" $NO_CACHE $BUILD_ARGS .
else
  echo "Building Docker image with multi-stage build (this will take several minutes)..."
  docker build -f Dockerfile -t "$TAG" $NO_CACHE $BUILD_ARGS .
fi

echo "Docker image built successfully: $TAG"
echo ""
echo "Usage examples:"
echo "  docker run --rm $TAG STAR --help"
echo "  docker run --rm -v \$(pwd)/data:/data $TAG STAR --genomeDir /data/genome --readFilesIn /data/reads.fastq"
