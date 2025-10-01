#!/usr/bin/env python3
"""
Analyze concordance between ZG (new) and GX (original) gene tags in BAM file.
"""

import sys
import pysam
from collections import defaultdict

def main():
    if len(sys.argv) != 2:
        print("Usage: analyze_concordance.py <bam_file>")
        sys.exit(1)

    bam_path = sys.argv[1]
    
    # Statistics
    total_reads = 0
    reads_with_gx = 0
    reads_with_zg = 0
    reads_with_both = 0
    perfect_matches = 0
    zg_superset = 0  # ZG contains all GX genes plus more
    gx_superset = 0  # GX contains all ZG genes plus more  
    partial_overlap = 0
    no_overlap = 0
    
    # Detailed comparison
    gene_comparisons = defaultdict(int)
    
    print(f"Analyzing concordance in: {bam_path}")
    print("=" * 60)
    
    try:
        with pysam.AlignmentFile(bam_path, 'rb') as bam:
            for rec in bam:
                total_reads += 1
                
                # Get tags
                gx_value = rec.get_tag('GX') if rec.has_tag('GX') else None
                zg_value = rec.get_tag('ZG') if rec.has_tag('ZG') else None
                
                # Parse gene sets
                gx_genes = set()
                zg_genes = set()
                
                if gx_value and gx_value != '-':
                    reads_with_gx += 1
                    gx_genes = set(gx_value.split(','))
                
                if zg_value and zg_value != '-':
                    reads_with_zg += 1
                    zg_genes = set(zg_value.split(','))
                
                if gx_genes and zg_genes:
                    reads_with_both += 1
                
                # Compare gene sets
                if gx_genes or zg_genes:
                    if gx_genes == zg_genes:
                        perfect_matches += 1
                        gene_comparisons['perfect_match'] += 1
                    elif zg_genes.issuperset(gx_genes):
                        zg_superset += 1
                        gene_comparisons['zg_superset'] += 1
                    elif gx_genes.issuperset(zg_genes):
                        gx_superset += 1
                        gene_comparisons['gx_superset'] += 1
                    elif gx_genes & zg_genes:  # Some overlap
                        partial_overlap += 1
                        gene_comparisons['partial_overlap'] += 1
                    else:  # No overlap
                        no_overlap += 1
                        gene_comparisons['no_overlap'] += 1
                
                # Show first 10 examples for debugging
                if total_reads <= 10 and (gx_genes or zg_genes):
                    print(f"Read {total_reads}: GX={gx_genes}, ZG={zg_genes}")

    except Exception as e:
        print(f"Error processing BAM file: {e}")
        sys.exit(1)

    # Print results
    print(f"\n=== Concordance Analysis Results ===")
    print(f"Total reads processed: {total_reads:,}")
    print(f"Reads with GX tags: {reads_with_gx:,} ({reads_with_gx/total_reads*100:.1f}%)")
    print(f"Reads with ZG tags: {reads_with_zg:,} ({reads_with_zg/total_reads*100:.1f}%)")
    print(f"Reads with both tags: {reads_with_both:,} ({reads_with_both/total_reads*100:.1f}%)")
    
    print(f"\n=== Gene Set Comparisons ===")
    annotated_reads = sum(gene_comparisons.values())
    if annotated_reads > 0:
        print(f"Total annotated reads: {annotated_reads:,}")
        print(f"Perfect matches (GX == ZG): {perfect_matches:,} ({perfect_matches/annotated_reads*100:.1f}%)")
        print(f"ZG superset (ZG ⊃ GX): {zg_superset:,} ({zg_superset/annotated_reads*100:.1f}%)")
        print(f"GX superset (GX ⊃ ZG): {gx_superset:,} ({gx_superset/annotated_reads*100:.1f}%)")
        print(f"Partial overlap: {partial_overlap:,} ({partial_overlap/annotated_reads*100:.1f}%)")
        print(f"No overlap: {no_overlap:,} ({no_overlap/annotated_reads*100:.1f}%)")
        
        # Overall concordance metric
        high_concordance = perfect_matches + zg_superset  # ZG should be superset due to GeneFull
        print(f"\nOverall concordance: {high_concordance:,}/{annotated_reads:,} ({high_concordance/annotated_reads*100:.1f}%)")
        
        if zg_superset > perfect_matches:
            print("✓ ZG tags capture more genes than GX (expected with GeneFull)")
        elif perfect_matches > zg_superset:
            print("✓ ZG and GX tags are highly concordant")
        else:
            print("⚠ Mixed results - investigate discrepancies")
    else:
        print("No annotated reads found")

if __name__ == "__main__":
    main()
