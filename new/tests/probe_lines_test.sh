#!/usr/bin/env bash
set -euo pipefail

# Probe-line integration test: validates that BAM contains records whose RNAME (col 3)
# does not start with 'chr' (probe references), and inspects unassigned and multimapper cases.

if [[ $# -lt 1 ]]; then
  echo "Usage: $(basename "$0") <Aligned.out.bam>" >&2
  exit 2
fi

BAM_PATH=$1

if ! command -v samtools >/dev/null 2>&1; then
  echo "samtools not found in PATH" >&2
  exit 1
fi

if [[ ! -f "$BAM_PATH" ]]; then
  echo "BAM not found: $BAM_PATH" >&2
  exit 1
fi

echo "BAM: $BAM_PATH"

# Count probe reference lines (RNAME not starting with 'chr' and not '*')
probe_count=$(samtools view "$BAM_PATH" | awk '($3 !~ /^chr/ && $3 != "*") {c++} END{print c+0}')
echo "Probe-reference records (RNAME !~ ^chr): $probe_count"

# Sample 5 probe lines and print RNAME, NH, ZG, GX
echo "Sample probe lines (RNAME, NH, ZG, GX):"
samtools view "$BAM_PATH" | awk '($3 !~ /^chr/ && $3 != "*") {print $3,$0}' | head -n 5 | \
  awk '{rname=$1; $1=""; line=$0; nh="NA"; zg="NA"; gx="NA"; split(line, f, /\t/); for(i=1;i<=length(f);i++){if(f[i]~"^NH:i:") nh=substr(f[i],6); if(f[i]~"^ZG:Z:") zg=substr(f[i],6); if(f[i]~"^GX:Z:") gx=substr(f[i],6);} printf("  %s\tNH=%s\tZG=%s\tGX=%s\n", rname, nh, zg, gx)}'

# Count multimappers among probe lines (NH>1)
probe_multi=$(samtools view "$BAM_PATH" | awk '($3 !~ /^chr/ && $3 != "*") {if(match($0,/NH:i:([0-9]+)/,m)&&m[1]>1) c++} END{print c+0}')
echo "Probe multimappers (NH>1): $probe_multi"

# Count unassigned among all lines (GX:Z:- or GN:Z:-)
unassigned=$(samtools view "$BAM_PATH" | awk '{if($0 ~ /GX:Z:-/ || $0 ~ /GN:Z:-/) c++} END{print c+0}')
echo "Unassigned (GX:Z:- or GN:Z:-): $unassigned"

# Basic assertions: require at least 1 probe line; do not fail if multimappers/unassigned are zero
if [[ "$probe_count" -lt 1 ]]; then
  echo "No probe-reference records found; check input BAM." >&2
  exit 2
fi

echo "Probe-line test passed."

