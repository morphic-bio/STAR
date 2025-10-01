#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: build_star.sh [options]

Build STAR from source using the Makefile in source/.

Options:
  --static                Build static binary variants (STAR, optionally STARlong)
  --long                  Also build long-reads variant (STARlong)
  --only-long             Build only the long-reads variant (STARlong)
  --debug                 Build debug variant (gdb or gdb-long); ignores --static
  --posix-shared          Build with POSIX shared memory support (POSIXSHARED target)
  --simd <flag>           Override SIMD flag passed to compiler (e.g. -msse4.2, -mavx, -mavx2)
  --cxx <path>            Path to C++ compiler (g++)
  --cxxflags-extra <str>  Extra CXX flags to append
  --ldflags-extra <str>   Extra LDFLAGS to append
  --jobs, -j <N>          Parallel jobs for make (default: detected core count)
  --clean                 Run a clean build (invokes "make CLEAN")
  --no-copy               Do not copy built binaries into bin/
  --prefix <dir>          Copy built binaries to this directory instead of bin/
  --target <name>         Build a specific Makefile target (advanced; bypasses other flags)
  -h, --help              Show this help

Examples:
  ./build_star.sh                         # Build default STAR (dynamic)
  ./build_star.sh --long                  # Build STAR and STARlong (dynamic)
  ./build_star.sh --static --long         # Build static STAR and STARlong
  ./build_star.sh --simd -msse4.2         # Build with SSE4.2 target
  ./build_star.sh --cxx /usr/bin/g++-12   # Use a specific compiler
EOF
}

msg() { echo "[build_star] $*"; }

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$REPO_ROOT/source"

if [[ ! -d "$SOURCE_DIR" ]]; then
  echo "Error: source directory not found at $SOURCE_DIR" >&2
  exit 1
fi

STATIC=0
BUILD_LONG=0
ONLY_LONG=0
DEBUG=0
POSIX_SHARED=0
CLEAN=0
COPY_TO_BIN=1
JOBS=""
CXX_BIN=""
SIMD_FLAG=""
PREFIX_DIR=""
LDFLAGS_EXTRA=""
CXXFLAGS_EXTRA=""
EXPLICIT_TARGET=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --static) STATIC=1 ;;
    --long) BUILD_LONG=1 ;;
    --only-long) ONLY_LONG=1; BUILD_LONG=1 ;;
    --debug) DEBUG=1 ;;
    --posix-shared) POSIX_SHARED=1 ;;
    --clean) CLEAN=1 ;;
    --no-copy) COPY_TO_BIN=0 ;;
    --jobs|-j) JOBS="${2:-}"; shift ;;
    --cxx) CXX_BIN="${2:-}"; shift ;;
    --simd) SIMD_FLAG="${2:-}"; shift ;;
    --prefix) PREFIX_DIR="${2:-}"; shift ;;
    --ldflags-extra) LDFLAGS_EXTRA="${2:-}"; shift ;;
    --cxxflags-extra) CXXFLAGS_EXTRA="${2:-}"; shift ;;
    --target) EXPLICIT_TARGET="${2:-}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
  shift
done

# Derive default job count
if [[ -z "$JOBS" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    JOBS="$(nproc)"
  elif [[ "$(uname -s)" == "Darwin" ]]; then
    JOBS="$(sysctl -n hw.ncpu)"
  else
    JOBS=1
  fi
fi

# Resolve bin subdir based on platform and static/dynamic
OS_NAME="$(uname -s)"
ARCH_NAME="$(uname -m)"
case "${OS_NAME}/${ARCH_NAME}" in
  Linux/x86_64) BIN_OS_DIR="Linux_x86_64" ;;
  Darwin/x86_64) BIN_OS_DIR="MacOSX_x86_64" ;;
  *) BIN_OS_DIR="${OS_NAME}_${ARCH_NAME}" ;;
esac
if (( STATIC )); then
  BIN_OS_DIR="${BIN_OS_DIR}_static"
fi

# Compose make variable overrides
MAKE_VARS=()
if [[ -n "$CXX_BIN" ]]; then
  MAKE_VARS+=("CXX=$CXX_BIN")
fi
if [[ -n "$SIMD_FLAG" ]]; then
  MAKE_VARS+=("CXXFLAGS_SIMD=$SIMD_FLAG")
fi
if [[ -n "$CXXFLAGS_EXTRA" ]]; then
  MAKE_VARS+=("CXXFLAGSextra=$CXXFLAGS_EXTRA")
fi
if [[ -n "$LDFLAGS_EXTRA" ]]; then
  MAKE_VARS+=("LDFLAGSextra=$LDFLAGS_EXTRA")
fi

# Decide targets to build
TARGETS=()
if [[ -n "$EXPLICIT_TARGET" ]]; then
  TARGETS+=("$EXPLICIT_TARGET")
else
  if (( DEBUG )); then
    if (( BUILD_LONG )); then
      TARGETS+=("gdb-long")
    else
      TARGETS+=("gdb")
    fi
  elif (( POSIX_SHARED )); then
    TARGETS+=("POSIXSHARED")
  else
    if (( STATIC )); then
      if (( ONLY_LONG )); then
        TARGETS+=("STARlongStatic")
      else
        TARGETS+=("STARstatic")
        if (( BUILD_LONG )); then TARGETS+=("STARlongStatic"); fi
      fi
    else
      if (( ONLY_LONG )); then
        TARGETS+=("STARlong")
      else
        TARGETS+=("STAR")
        if (( BUILD_LONG )); then TARGETS+=("STARlong"); fi
      fi
    fi
  fi
fi

if (( CLEAN )); then
  msg "Cleaning previous build artifacts..."
  ( cd "$SOURCE_DIR" && make CLEAN )
fi

msg "Building targets: ${TARGETS[*]} (jobs: $JOBS)"

for tgt in "${TARGETS[@]}"; do
  ( cd "$SOURCE_DIR" && make -j"$JOBS" "$tgt" "${MAKE_VARS[@]}" )

  # Determine produced binary name
  case "$tgt" in
    *long*) BIN_NAME="STARlong" ;;
    *) BIN_NAME="STAR" ;;
  esac

  BIN_PATH="$SOURCE_DIR/$BIN_NAME"
  if [[ ! -x "$BIN_PATH" ]]; then
    echo "Error: expected binary not found: $BIN_PATH" >&2
    exit 1
  fi

  if (( COPY_TO_BIN )); then
    if [[ -n "$PREFIX_DIR" ]]; then
      DEST_DIR="$PREFIX_DIR"
    else
      DEST_DIR="$REPO_ROOT/bin/$BIN_OS_DIR"
    fi
    mkdir -p "$DEST_DIR"
    install -m 0755 "$BIN_PATH" "$DEST_DIR/$BIN_NAME"
    msg "Installed $BIN_NAME -> $DEST_DIR/$BIN_NAME"
  else
    msg "Built $BIN_PATH (not copied)"
  fi
done

msg "Done."


