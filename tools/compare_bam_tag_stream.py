#!/usr/bin/env python3
"""
Compare BAM CB/UB tags with binary tag stream for validation.

This script performs a deterministic, streaming comparison that walks the BAM 
and the serialized tag stream in lockstep and validates CB/UB presence and 
values per record.
"""

import sys
import argparse
import struct
import subprocess
import math
from typing import List, Tuple, Optional

def read_whitelist(path: str) -> List[str]:
    """Read CB whitelist from file, one per line."""
    whitelist = []
    try:
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    whitelist.append(line)
    except Exception as e:
        print(f"ERROR: Cannot read whitelist file {path}: {e}", file=sys.stderr)
        sys.exit(1)
    return whitelist

def read_binary_header(path: str) -> Tuple[int, int, int, int]:
    """Read binary header: statusBits, cbBits, umiBits, recordCount."""
    try:
        with open(path, 'rb') as f:
            header_data = f.read(32)
            if len(header_data) < 32:
                print(f"ERROR: Binary header truncated (got {len(header_data)} bytes, expected 32)", file=sys.stderr)
                sys.exit(1)
            
            # Unpack 4x uint64_t in little-endian format
            status_bits, cb_bits, umi_bits, record_count = struct.unpack('<QQQQ', header_data)
            return status_bits, cb_bits, umi_bits, record_count
    except Exception as e:
        print(f"ERROR: Cannot read binary header from {path}: {e}", file=sys.stderr)
        sys.exit(1)

def decode_umi(packed_umi: int, umi_length: int) -> str:
    """Decode packed UMI to string. MSB -> first base."""
    if packed_umi == 0:
        return "-"
    
    bases = ['A', 'C', 'G', 'T']
    umi = []
    
    # Read bits from high order end (MSB to LSB) to match STAR's packing order
    for i in range(umi_length):
        base_idx = (packed_umi >> (2 * (umi_length - 1 - i))) & 3
        umi.append(bases[base_idx])
    
    return ''.join(umi)

def parse_bam_record(sam_line: str) -> Tuple[str, str]:
    """Extract CB and UB from SAM line. Returns (cb, ub) or ("-", "-") if missing."""
    fields = sam_line.strip().split('\t')
    if len(fields) < 12:
        return "-", "-"
    
    cb = "-"
    ub = "-"
    
    # Parse optional fields (starting from field 12)
    for field in fields[11:]:  # 0-indexed, so field 11 is the 12th field
        if field.startswith('CB:Z:'):
            cb = field[5:]  # Remove 'CB:Z:' prefix
        elif field.startswith('UB:Z:'):
            ub = field[5:]  # Remove 'UB:Z:' prefix
    
    return cb, ub

def stream_bam_records(bam_path: str):
    """Stream BAM records via samtools view."""
    try:
        proc = subprocess.Popen(['samtools', 'view', bam_path], 
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE,
                               text=True)
        
        for line in proc.stdout:
            yield line
        
        proc.wait()
        if proc.returncode != 0:
            stderr_output = proc.stderr.read()
            print(f"ERROR: samtools view failed: {stderr_output}", file=sys.stderr)
            sys.exit(1)
            
    except Exception as e:
        print(f"ERROR: Cannot stream BAM file {bam_path}: {e}", file=sys.stderr)
        sys.exit(1)

def read_binary_records(path: str, record_count: int, cb_bits: int, umi_bits: int):
    """Generator that yields (status, cb_idx, umi_packed) for each record."""
    record_bits = 1 + cb_bits + umi_bits  # status + CB + UMI
    record_bytes = math.ceil(record_bits / 8)
    
    try:
        with open(path, 'rb') as f:
            # Skip header (32 bytes)
            f.seek(32)
            
            for i in range(record_count):
                record_data = f.read(record_bytes)
                if len(record_data) < record_bytes:
                    print(f"ERROR: Binary record {i} truncated (got {len(record_data)} bytes, expected {record_bytes})", file=sys.stderr)
                    sys.exit(1)
                
                # Convert bytes to integer (little-endian)
                value = 0
                for j, byte in enumerate(record_data):
                    value |= byte << (j * 8)
                
                # Extract fields (LSB-first)
                status = value & ((1 << 1) - 1)
                value >>= 1
                
                cb_idx = value & ((1 << cb_bits) - 1)
                value >>= cb_bits
                
                umi_packed = value & ((1 << umi_bits) - 1)
                
                yield status, cb_idx, umi_packed
                
    except Exception as e:
        print(f"ERROR: Cannot read binary records from {path}: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Compare BAM CB/UB tags with binary tag stream')
    parser.add_argument('--bam', required=True, help='BAM file path')
    parser.add_argument('--tags', required=True, help='Binary tag stream path')
    parser.add_argument('--whitelist', required=True, help='CB whitelist file path')
    parser.add_argument('--limit', type=int, help='Limit comparison to first N records (for testing)')
    
    args = parser.parse_args()
    
    # Read whitelist
    whitelist = read_whitelist(args.whitelist)
    print(f"Loaded {len(whitelist)} CB entries from whitelist", file=sys.stderr)
    
    # Read binary header
    status_bits, cb_bits, umi_bits, record_count = read_binary_header(args.tags)
    print(f"Binary header: statusBits={status_bits}, cbBits={cb_bits}, umiBits={umi_bits}, recordCount={record_count}", file=sys.stderr)
    
    # Validate header
    if status_bits != 1:
        print(f"ERROR: Expected statusBits=1, got {status_bits}", file=sys.stderr)
        sys.exit(1)
    
    if umi_bits % 2 != 0:
        print(f"ERROR: umiBits must be even (got {umi_bits})", file=sys.stderr)
        sys.exit(1)
    
    umi_length = umi_bits // 2
    print(f"UMI length: {umi_length} bases", file=sys.stderr)
    
    # Apply limit if specified
    if args.limit:
        record_count = min(record_count, args.limit)
        print(f"Limiting comparison to first {record_count} records", file=sys.stderr)
    
    # Compare records
    matches = 0
    errors = 0
    
    bam_stream = stream_bam_records(args.bam)
    binary_stream = read_binary_records(args.tags, record_count, cb_bits, umi_bits)
    
    for record_idx, (binary_record, bam_line) in enumerate(zip(binary_stream, bam_stream)):
        if record_idx >= record_count:
            break
            
        status, cb_idx, umi_packed = binary_record
        bam_cb, bam_ub = parse_bam_record(bam_line)
        
        # Comparison logic based on status
        if status == 0:
            # Status=0: UB must be missing/empty, CB ignored
            if bam_ub not in ["-", ""]:
                print(f"ERROR: Record {record_idx}: Binary status=0 but BAM has UB='{bam_ub}'", file=sys.stderr)
                errors += 1
                if errors >= 10:  # Limit error output
                    print("... (stopping after 10 errors)", file=sys.stderr)
                    break
            else:
                matches += 1
        else:
            # Status=1: Both CB and UB should be present and match
            if cb_idx == 0 or cb_idx > len(whitelist):
                print(f"ERROR: Record {record_idx}: Invalid CB index {cb_idx} (whitelist size: {len(whitelist)})", file=sys.stderr)
                errors += 1
                if errors >= 10:
                    print("... (stopping after 10 errors)", file=sys.stderr)
                    break
                continue
            
            expected_cb = whitelist[cb_idx - 1]  # -1 because CB indices are 1-based
            expected_ub = decode_umi(umi_packed, umi_length)
            
            if bam_cb != expected_cb:
                print(f"ERROR: Record {record_idx}: CB mismatch (BAM: '{bam_cb}', Expected: '{expected_cb}')", file=sys.stderr)
                errors += 1
            elif bam_ub != expected_ub:
                print(f"ERROR: Record {record_idx}: UB mismatch (BAM: '{bam_ub}', Expected: '{expected_ub}')", file=sys.stderr)
                errors += 1
            else:
                matches += 1
            
            if errors >= 10:
                print("... (stopping after 10 errors)", file=sys.stderr)
                break
    
    # Check for record count mismatch
    remaining_bam = list(bam_stream)
    if len(remaining_bam) > 0:
        print(f"ERROR: BAM has {len(remaining_bam)} extra records after {record_count} binary records", file=sys.stderr)
        errors += 1
    
    # Summary
    total_checked = matches + errors
    if errors == 0:
        print(f"SUCCESS: All {matches} records match correctly")
        sys.exit(0)
    else:
        print(f"FAILURE: {errors} mismatches out of {total_checked} records checked", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()