#!/usr/bin/env bash

# Plan: Fix emit_test.sh CB/UB comparison using a single-pass, ordered, record-by-record check

# Objective
# Replace the brittle TSV-based joins and ad‑hoc decoders in emit_test.sh with a
# deterministic, streaming comparison that walks the BAM and the serialized tag
# stream in lockstep and validates CB/UB presence and values per record.

# High‑Level Approach
# - Add a Python helper: tools/compare_bam_tag_stream.py
#   • Reads 32‑byte binary header (LE) => statusBits, cbBits, umiBits, recordCount
#   • Computes bytes/record and iterates recordCount times, decoding per‑record
#     status, CB index (1‑based, 0 sentinel), and packed UMI
#   • Streams BAM alignments via `samtools view` and extracts CB:Z: / UB:Z:
#   • For each record i:
#       - If status==0: BAM UB must be missing/empty; CB ignored
#       - If status==1: CB index must be in whitelist; compare CB string, decode
#         UMI to string and compare to UB; both must match
#   • Verifies BAM record count == header.recordCount and exits non‑zero on any
#     mismatch or truncation
# - Update emit_test.sh to replace direct_bam_binary_comparison() body with a
#   simple call to the Python helper.

# Preconditions
# - samtools available on PATH (used to stream BAM as SAM)
# - python3 available
# - A CB whitelist path (existing WHITELIST variable in emit_test.sh)
# - The binary tag stream file produced by STAR (Aligned.out.cb_ub.bin)

# Implementation Steps

# 1) Add helper script tools/compare_bam_tag_stream.py
#    Content/specification:
#    - CLI args:
#        --bam <path>         BAM path with CB/UB injected (unsorted is fine)
#        --tags <path>        Binary tag stream path
#        --whitelist <path>   Text file of CB whitelist (1 per line)
#    - Header parsing:
#        * Read 32 bytes, LE unpack to Q Q Q Q => statusBits, cbBits, umiBits, recordCount
#        * Validate statusBits==1; compute record_bits = 1 + cbBits + umiBits;
#          record_bytes = ceil(record_bits/8)
#        * Derive umi_length = umiBits/2; enforce umiBits even
#    - Record decoding:
#        * For each i in [0, recordCount): read record_bytes, cast to little‑endian int
#        * Extract fields (LSB‑first): status = value & ((1<<1)-1); value >>=1
#                                      cbIdx = value & ((1<<cbBits)-1); value >>=cbBits
#                                      umiPacked = value & ((1<<umiBits)-1)
#    - UMI unpack:
#        * Convert packed 2‑bit bases to string of umi_length; correct orientation
#          is MSB→first base; implement by filling from end while consuming LSBs
#    - BAM streaming:
#        * Spawn `samtools view <bam>` with text pipe; skip header lines if present
#        * For each record, parse optional fields from column 12 onward; extract CB:Z: and UB:Z:
#    - Comparison rules:
#        * status==0: UB must be absent/empty/"-"; CB ignored. Mismatch => error
#        * status==1: cbIdx>0 and <=len(whitelist); expected_cb = whitelist[cbIdx-1]
#                     UB string must equal decoded UMI; CB string must equal expected_cb
#    - End conditions:
#        * If BAM ends before recordCount => error; if BAM has extra records after
#          recordCount => error
#        * On any mismatch, print diagnostic with record index and values and exit(1)
#        * On success, print a single summary line and exit(0)
#    - Performance notes:
#        * Process line‑by‑line; avoid buffering entire files; suitable for large runs
#
# 2) Make the helper executable and ensure it’s committed
#    - chmod +x tools/compare_bam_tag_stream.py
#    - Add to repo (handled by commit step outside this plan)
#
# 3) Replace direct_bam_binary_comparison() in emit_test.sh
#    - Locate the function definition:
#        direct_bam_binary_comparison() { ... }
#    - Replace its body with a wrapper that calls the Python helper:
#        direct_bam_binary_comparison() {
#          local bam_file="$1"; local binary_file="$2"; local description="$3"
#          echo "Direct comparison: $description"; echo "  BAM: $bam_file"; echo "  Binary: $binary_file"
#          local compare_script="$(dirname "$0")/tools/compare_bam_tag_stream.py"
#          if [[ ! -x "$compare_script" ]]; then echo "ERROR: $compare_script not executable"; return 1; fi
#          if python3 "$compare_script" --bam "$bam_file" --tags "$binary_file" --whitelist "$WHITELIST"; then
#            echo "  ✓ BAM and tag stream are consistent"; return 0; else
#            echo "  ✗ BAM and tag stream comparison failed"; return 1; fi
#        }
#    - Remove any now‑unused helpers in emit_test.sh such as:
#        • extract_bam_cb_ub_with_record_index
#        • extract_table_cb_ub
#        • The TSV join logic under the previous comparison implementation
#    - Keep validate_tag_table as a sanity check if desired, or simplify it to only
#      call the Python helper in a dry‑run mode (optional).
#
# 4) Update decoder orientation (optional cleanup)
#    - If tools/decode_tag_binary is still used elsewhere, ensure its UMI unpacking
#      matches STAR’s encoding (MSB→leftmost base). The new Python helper already
#      implements the correct orientation; consider deprecating the C++ decoder or
#      fixing it in place to avoid future confusion.
#
# 5) Parameterization
#    - UMI length: derive from umiBits/2 in the Python helper; do not hard‑code 12
#    - Whitelist: re‑use existing WHITELIST var in emit_test.sh; pass through
#    - Endianness: the helper assumes little‑endian header and LSB‑first packing,
#      which matches BAMTagBinaryWriter
#
# 6) Failure Modes & Diagnostics
#    - Truncated header (len<32) => fail with message
#    - statusBits!=1 or odd umiBits => fail with message
#    - Truncated record (short read) => fail with message and record index
#    - BAM exhausted early or has extra records => fail with counts and indices
#    - CB index out of range => fail with index value
#    - CB/UB mismatches => fail with record index and show BAM vs expected values
#
# 7) CI/Test Integration
#    - Ensure emit_test.sh exits non‑zero if direct comparison fails
#    - Optionally add a fast mode that samples first N records for smoke tests
#    - Document that emit_test.sh requires samtools and python3
#
# 8) Manual Validation
#    - Run Test 2 (both modes) and confirm:
#        • validate_tag_table passes
#        • direct_bam_binary_comparison passes
#        • Solo matrices match baseline
#    - Spot‑check first 10 decoded records by printing from the Python helper
#      (add a --limit option temporarily if useful)
#
# 9) Clean‑up
#    - Remove obsolete temporary files and helper functions from emit_test.sh
#    - Keep tools/compare_bam_tag_stream.py as the single source of truth for
#      decoding/validation in test workflows

# End of plan

