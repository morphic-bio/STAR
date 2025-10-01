#!/usr/bin/env python3
"""
Final validation script for ZG/ZX BAM tags implementation.
Validates that ZG (gene set) and ZX (overlap status) tags are correctly formatted and contain expected values.
"""

import sys
import os
import pysam

def main():
    if len(sys.argv) != 3:
        print("Usage: validate_zg_zx.py <bam_file> <allowed_genes_file>")
        print("Validates ZG and ZX tags in BAM file against allowed gene list")
        sys.exit(1)

    bam_path = sys.argv[1]
    allowed_genes_path = sys.argv[2]

    # Load allowed genes
    allowed = set()
    if os.path.exists(allowed_genes_path):
        with open(allowed_genes_path, 'r') as f:
            allowed = {line.strip() for line in f if line.strip()}

    print(f"Validating BAM: {bam_path}")
    print(f"Allowed genes: {allowed if allowed else 'None (empty gene list)'}")

    # Validation counters
    total_reads = 0
    reads_with_zg = 0
    reads_with_zx = 0
    reads_with_genes = 0
    reads_with_overlap = 0
    
    valid_zx_values = {'none', 'exonic', 'intronic', 'intergenic', 'spanning'}

    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            for rec in bam:
                total_reads += 1
                
                # Validate ZG tag
                if not rec.has_tag('ZG'):
                    raise SystemExit(f"ERROR: Missing ZG tag for read {rec.query_name}")
                
                reads_with_zg += 1
                zg_value = rec.get_tag('ZG')
                
                if zg_value != '-':  # If not empty gene set
                    reads_with_genes += 1
                    zg_genes = set(zg_value.split(','))
                    
                    # Validate all genes are in allowed set (if allowed set is not empty)
                    if allowed and not zg_genes.issubset(allowed):
                        unexpected = zg_genes - allowed
                        raise SystemExit(f"ERROR: Unexpected genes in {rec.query_name}: {unexpected}")
                
                # Validate ZX tag
                if not rec.has_tag('ZX'):
                    raise SystemExit(f"ERROR: Missing ZX tag for read {rec.query_name}")
                
                reads_with_zx += 1
                zx_value = rec.get_tag('ZX')
                
                if zx_value not in valid_zx_values:
                    raise SystemExit(f"ERROR: Invalid ZX value '{zx_value}' for read {rec.query_name}. Valid values: {valid_zx_values}")
                
                if zx_value != 'none':
                    reads_with_overlap += 1

    except FileNotFoundError:
        raise SystemExit(f"ERROR: BAM file not found: {bam_path}")
    except Exception as e:
        raise SystemExit(f"ERROR: Failed to process BAM file: {e}")

    # Print validation results
    print(f"\n=== ZG/ZX Tag Validation Results ===")
    print(f"Total reads processed: {total_reads}")
    print(f"Reads with ZG tags: {reads_with_zg}")
    print(f"Reads with ZX tags: {reads_with_zx}")
    print(f"Reads with gene assignments: {reads_with_genes}")
    print(f"Reads with overlap annotations: {reads_with_overlap}")
    
    if total_reads == 0:
        raise SystemExit("ERROR: No reads found in BAM file")
    
    if reads_with_zg != total_reads:
        raise SystemExit(f"ERROR: Not all reads have ZG tags ({reads_with_zg}/{total_reads})")
    
    if reads_with_zx != total_reads:
        raise SystemExit(f"ERROR: Not all reads have ZX tags ({reads_with_zx}/{total_reads})")

    print("SUCCESS: All ZG/ZX tags are valid and properly formatted!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
