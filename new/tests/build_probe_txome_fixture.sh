#!/usr/bin/env bash
set -euo pipefail

# Build a tiny STAR transcriptome genomeDir from a subset of a probe-only FASTA/GTF.
# Produces a directory suitable for STAR's Transcriptome loader (geneInfo.tab, exonInfo.tab, etc.).
#
# Usage:
#   build_probe_txome_fixture.sh --fasta /path/to/probes_only.fa \
#                                --gtf /path/to/probes_only.gtf \
#                                --out /tmp/probe_txome \
#                                [--contigs name1,name2,...] [--n 2]
#
# After completion:
#   export STAR_TXOME_DIR=/tmp/probe_txome
#   cd source && make test

FASTA=""
GTF=""
OUTDIR=""
CONTIGS=""
NUM=2
THREADS=${THREADS:-1}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) FASTA="$2"; shift 2;;
    --gtf) GTF="$2"; shift 2;;
    --out) OUTDIR="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --n) NUM="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help)
      sed -n '1,40p' "$0"; exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 2;;
  esac
done

if [[ -z "$FASTA" || -z "$GTF" || -z "$OUTDIR" ]]; then
  echo "Usage: $0 --fasta <fa> --gtf <gtf> --out <dir> [--contigs name1,name2] [--n 2]" >&2
  exit 2
fi

if [[ ! -f "$FASTA" ]]; then echo "FASTA not found: $FASTA" >&2; exit 1; fi
if [[ ! -f "$GTF" ]]; then echo "GTF not found: $GTF" >&2; exit 1; fi

mkdir -p "$OUTDIR"
WORKDIR=$(mktemp -d "${OUTDIR%/}/build.XXXXXX")
trap 'rm -rf "$WORKDIR"' EXIT

# Resolve STAR binary: prefer repo-local build, fall back to PATH
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
STAR_BIN="$REPO_ROOT/source/STAR"
if [[ ! -x "$STAR_BIN" ]]; then
  if command -v STAR >/dev/null 2>&1; then
    STAR_BIN=$(command -v STAR)
  else
    echo "STAR binary not found at $REPO_ROOT/source/STAR and not in PATH. Build it first (cd source && make STAR)." >&2
    exit 1
  fi
fi

# Decide which contigs to keep
if [[ -z "$CONTIGS" ]]; then
  # Take top N contig names that exist in both FASTA and GTF
  fasta_names=$(awk '/^>/{split($1,h,">"); print h[2]}' "$FASTA" | head -n $NUM)
  gtf_names=$(awk '{print $1}' "$GTF" | sort -u)
  # Join sets, keep those present in GTF too
  CONTIGS=$(awk 'NR==FNR{a[$0]=1; next} ($0 in a){print}' <(echo "$fasta_names") <(echo "$gtf_names") | head -n $NUM | paste -sd, -)
fi

if [[ -z "$CONTIGS" ]]; then
  echo "Could not determine contigs to extract. Provide --contigs name1,name2" >&2
  exit 1
fi

echo "Using contigs: $CONTIGS"

MINI_FA="$WORKDIR/mini.fa"
MINI_GTF="$WORKDIR/mini.gtf"

# Extract selected contigs from FASTA
awk -v list="$CONTIGS" '
  BEGIN{n=split(list,a,","); for(i=1;i<=n;i++){sel[a[i]]=1}}
  /^>/{split(substr($0,2),t,/\s+/); name=t[1]; printing=(name in sel)}
  { if (printing) print }
' "$FASTA" > "$MINI_FA"

# Filter GTF rows by selected contigs
awk -v list="$CONTIGS" '
  BEGIN{n=split(list,a,","); for(i=1;i<=n;i++){sel[a[i]]=1}}
  sel[$1]==1 { print }
' "$GTF" > "$MINI_GTF"

echo "Building STAR genomeDir at $OUTDIR"
"$STAR_BIN" --runMode genomeGenerate \
  --genomeDir "$OUTDIR" \
  --genomeFastaFiles "$MINI_FA" \
  --sjdbGTFfile "$MINI_GTF" \
  --sjdbOverhang 100 \
  --runThreadN "$THREADS" \
  --genomeSAindexNbases 6 \
  --limitGenomeGenerateRAM 1200000000

echo "Fixture built. Export and run tests:"
echo "  export STAR_TXOME_DIR=$OUTDIR"
echo "  cd $REPO_ROOT/source && make test"

