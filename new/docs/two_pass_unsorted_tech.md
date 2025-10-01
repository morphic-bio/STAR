# STARsolo Two-Pass Unsorted CB/UB Injection – Technical Notes

This document describes how the new unsorted CB/UB workflow is implemented and highlights code-level changes to STAR.

## Architecture Summary

```
Pass 1 (mapping threads)
  ├── ReadAlign::alignBAM() produces BAM records
  ├── If --soloAddTagsToUnsorted yes:
  │     └── BAMoutputSoloTmp::unsortedOneAlign()
  │            writes [uint32 size][payload][uint64 trailer] to a shared tmp file
  ├── If --soloWriteTagTable enabled:
  │     └── Standard BAMoutput (no temp file needed)
  └── Solo processes reads (unchanged)

Pass 2 (single-threaded post-processing)
  ├── If --soloAddTagsToUnsorted yes:
  │     └── BAMunsortedAddSoloTags()
  │           ▸ sequentially replays tmp file
  │           ▸ appends corrected CB/UB via SoloFeature::addBAMtags()
  │           ▸ writes final Aligned.out.bam (unsorted) and deletes tmp file
  └── If --soloWriteTagTable enabled:
        └── SoloFeature::writeTagTableIfRequested()
              ▸ exports CB/UB assignments to sidecar TSV file
              ▸ maintains strict alignment with BAM record order
```

## Key Modules

### 1. Parameter Handling (`source/Parameters*.{h,cpp}`)
- Added `ParametersSolo::{addTagsToUnsortedStr, addTagsToUnsorted}` with yes/no validation.
- Added `ParametersSolo::{writeTagTableStr, writeTagTableEnabled, writeTagTablePath}` for tag table export.
- Hooked the new CLI flags in `Parameters.cpp` and `parametersDefault`.
- After `pSolo.initialize(this)` runs, the unsorted BAM code path decides between:
  - opening the traditional BGZF stream (`outBAMfileUnsorted`), or
  - opening a shared `ofstream` for the solo tmp (`outBAMfileUnsortedSoloTmp`).
- Tag table mode bypasses the temp BAM entirely unless `--soloAddTagsToUnsorted` is also enabled.

### 2. Pass-1 Capture (`source/BAMoutputSoloTmp.{h,cpp}`)
- New lightweight writer that accepts a BAM record plus `iReadAll`, buffering output with mutex protection.
- Guardrails ensure any record exceeding `BAMoutput_oneAlignMaxBytes` aborts cleanly.
- Invoked from `ReadAlign_outputAlignments.cpp` and the chimeric writer; execution order now checks `outBAMsoloTmp` **before** the legacy unsorted path to avoid null dereferences.

### 3. Pass-2 Replay (`source/BAMunsortedAddSoloTags.{h,cpp}`)
- Streams the tmp file record-by-record. For each entry:
  1. Validate payload size.
  2. Append the trailer (containing `iReadAll`) to the buffer.
  3. Call `SoloFeature::addBAMtags` to append CB/UB tags.
  4. Write the resulting record into a freshly opened BGZF unsorted BAM.
- On success the temporary file is deleted; warnings are logged if removal fails.

### 4. STAR Orchestration (`source/STAR.cpp`)
- After `soloMain.processAndOutput()`, if the flag is active, call `BAMunsortedAddSoloTags` with the tmp path and desired BAM path.
- Logging announces pass-2 start/finish for easier monitoring.
- Cleanup closes the shared tmp stream prior to replay.

### 5. ReadAlign Changes
- `ReadAlign` now carries an `outBAMsoloTmp` pointer.
- All unsorted write sites (mapped, unmapped, chimeric) branch to the solo tmp writer when present.

### 5. Tag Table Export (`source/SoloFeature_writeTagTable.cpp`)
- New method `SoloFeature::writeTagTableIfRequested()` exports CB/UB assignments to a sidecar TSV file.
- Called after `collapseUMIall()` completes but before `outputResults()` in `SoloFeature::processRecords()`.
- Only runs when `writeTagTableEnabled` is true and `featureType == samAttrFeature` (usually Gene).
- Outputs tab-separated columns: `bam_record_index`, `iReadAll`, `mate`, `align_idx`, `qname`, `CB`, `UB`, `status`.
- Uses existing `readInfo` structure and `cbWLstr`/`convertNuclInt64toString` for CB/UMI string conversion.
- Memory-efficient: clears `readNames` vector after export to reduce peak RAM usage.

## Testing Additions
- `testing/scripts/10_run_sorted_baseline.sh`: runs the canonical sorted pipeline for comparison.
- `testing/scripts/20_run_unsorted_two_pass.sh`: executes the new two-pass mode.
- `testing/scripts/30_compare_bam_tags.py`, `31_compare_alignment_qc.py`, `32_compare_solo_matrices.sh`, `33_compare_logs.sh`, `34_compare_read_distributions.py`: confirm behavioural parity between runs.
- `testing/scripts/22_verify_tmp_trailers.py`: optional tmp-file integrity check (enabled with `KEEP_SOLO_TMP=1`).

## Data Layout of Solo Tmp File
Each entry consists of:
1. `uint32 size_with_len` – number of bytes in the BAM payload (excluding size field).
2. `size_with_len` bytes – raw BAM record as emitted by `ReadAlign::alignBAM`.
3. `uint64 trailer` – `iReadAll` stored in the upper 32 bits (lower 32 bits zeroed).

`SoloFeature::addBAMtags` expects this trailer immediately following the BAM record, so pass 2 `memcpy`s it back into place prior to tag injection.

## Buffer Safety
- `BAMoutputSoloTmp::unsortedOneAlign` checks `bamSize <= BAMoutput_oneAlignMaxBytes` before writing.
- `BAMunsortedAddSoloTags` re-validates `size_with_len` and the final size returned by `addBAMtags` against `BAM_ATTR_MaxSize`.
- Any violation aborts with a descriptive error, avoiding silent corruption.

## Notes on Concurrency
- The shared tmp stream is opened once in `Parameters.cpp`; each chunk instantiates its own writer but shares the underlying stream, protected by `g_threadChunks.mutexOutSAM` during writes.
- Flushing occurs on chunk completion and at shutdown; pass 2 runs after mapping threads finish, ensuring exclusive read access.

## Compatibility
- Legacy behaviour (unsorted without CB/UB) is unchanged.
- Sorted BAM output and all Solo features continue to work as before; the new flag only affects the unsorted path when CB/UB attributes are requested.

## Future Enhancements (Ideas)
- Optional BGZF compression of the tmp container to reduce disk usage.
- Multi-threaded pass 2 (sharded tmp files) for extreme datasets.
- Resume support by retaining the tmp file if pass 2 fails midway.

For end-user guidance, consult `docs/two_pass_unsorted_usage.md`.
