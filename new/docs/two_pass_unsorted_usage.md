# STARsolo Two-Pass Unsorted CB/UB Tag Injection – Usage Guide

This guide explains how to enable the new two-pass workflow that appends corrected CB/UB tags to **unsorted** BAM output without running the coordinate sorter.

## Prerequisites
- STAR genome index (e.g. `/storage/scRNAseq_output/indices-98-32/star`).
- FASTQ pairs (optionally lane-split) accessible to the run.
- 10x-compatible whitelist (e.g. `/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt`).
- Sufficient disk space for the temporary solo container (same order of magnitude as the final BAM).

## New/Updated Command-Line Options

| Flag | Default | Description |
|------|---------|-------------|
| `--soloAddTagsToUnsorted {yes|no}` | `no` | When `yes` and `--outSAMtype BAM Unsorted` is requested, STAR runs the two-pass workflow, generating a temporary file during pass 1 and producing an unsorted BAM with corrected CB/UB tags during pass 2. |
| `--soloWriteTagTable {None|Default|<path>}` | `None` | Export corrected CB/UB assignments to a sidecar table instead of or in addition to patching BAM. `None` disables (default), `Default` writes to `<outFileNamePrefix>Aligned.out.cb_ub.tsv`, custom path writes to specified file. |
| `--outSAMtype` | unchanged | Must include `BAM Unsorted`. The flag may still include `SortedByCoordinate`; the two-pass path only affects the `Unsorted` output. |

All other Solo options (barcode/UMI coordinates, features, filters) behave exactly as upstream STAR.

## Minimal Example

### Two-Pass BAM Tag Injection
```bash
STAR \
  --runThreadN 32 \
  --genomeDir /storage/scRNAseq_output/indices-98-32/star \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix /tmp/run1/ \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
  --outSAMtype BAM Unsorted \
  --soloAddTagsToUnsorted yes \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
```

### Tag Table Export (Alternative/Additional)
```bash
STAR \
  --runThreadN 32 \
  --genomeDir /storage/scRNAseq_output/indices-98-32/star \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix /tmp/run1/ \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt \
  --outSAMtype BAM Unsorted \
  --soloWriteTagTable Default \
  --outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM
```

### Outputs

#### Two-Pass Mode (`--soloAddTagsToUnsorted yes`)
- `Aligned.out.bam`: unsorted BAM containing corrected `CB:Z:` and `UB:Z:` tags.
- `Solo.out/`: standard Solo count matrices (unchanged format).
- Temporary file `Aligned.out.unsorted.solo.tmp` is created during pass 1 and removed after pass 2 succeeds. Set `KEEP_SOLO_TMP=1` in the environment to retain it for debugging.

#### Tag Table Mode (`--soloWriteTagTable Default`)
- `Aligned.out.bam`: standard unsorted BAM without CB/UB tags (unless `--soloAddTagsToUnsorted` is also enabled).
- `Aligned.out.cb_ub.tsv`: tab-separated sidecar table with corrected CB/UB assignments aligned with BAM records.
- `Solo.out/`: standard Solo count matrices (unchanged format).

#### Tag Table Format
The sidecar table contains the following columns:
1. `bam_record_index`: 0-based index incremented for every alignment record across mates and multimappers
2. `iReadAll`: 0-based STAR read counter (identical for both mates of a read)
3. `mate`: 0 for R1, 1 for R2
4. `align_idx`: 0-based index among multimappers for this read
5. `qname`: original read name
6. `CB`: corrected barcode string or `-` if unavailable
7. `UB`: corrected UMI string or `-` if unavailable
8. `status`: `OK`, `NO_CB`, `NO_UMI`, or `NO_CB_UMI`

## Workflow Overview
1. **Pass 1 – Capture**: Post-alignment BAC records are streamed into a temporary container via `BAMoutputSoloTmp`. Each entry stores the BAM record and the `iReadAll` identifier required for Solo correction.
2. **Solo Processing**: Solo aggregation/deduplication fills `readInfo[iread]` exactly as in the sorted pipeline.
3. **Pass 2 – Injection**: `BAMunsortedAddSoloTags` replays the temp file, calls `SoloFeature::addBAMtags`, and writes the final unsorted BAM. The temp file is deleted on success.

## Validation Checklist
Run the test scripts after producing both the sorted baseline and the unsorted run:
1. `bash testing/scripts/30_compare_bam_tags.py`
2. `bash testing/scripts/31_compare_alignment_qc.py`
3. `bash testing/scripts/32_compare_solo_matrices.sh`
4. `bash testing/scripts/33_compare_logs.sh`
5. `bash testing/scripts/34_compare_read_distributions.py`

For a turnkey end-to-end validation, execute `bash testing/scripts/90_run_all.sh`.

## Troubleshooting
- **Missing CB/UB tags in BAM**: ensure `--soloAddTagsToUnsorted yes` and that CB/UB are requested in `--outSAMattributes`.
- **Missing tag table**: ensure `--soloWriteTagTable` is not `None` and that `--soloFeatures` includes Gene/GeneFull.
- **Empty tag table**: check that CB/UB attributes are requested in `--outSAMattributes` to trigger Solo processing.
- **Temp file persists**: if STAR exits prematurely, inspect the temp file with `testing/scripts/22_verify_tmp_trailers.py`; delete it manually once done.
- **Large-record errors**: guardrails will abort if BAM records exceed the static buffer sizes. Adjust `BAMoutput_oneAlignMaxBytes` (in `IncludeDefine.h`) if necessary.

For deeper technical details, see `docs/two_pass_unsorted_tech.md`.
