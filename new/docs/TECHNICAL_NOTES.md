# Technical Notes – STARsolo Unsorted BAM & Tag Table Extensions

## Scope
This document captures the code changes introduced in the fork relative to upstream STAR 2.7.11b, focusing on unsorted BAM CB/UB injection and tag-table export. File paths refer to the repository root.

## Parameter & Configuration Flow
- `source/Parameters.cpp`: postpones unsorted BAM writer selection until after `pSolo.initialize()`, storing the decision in `outBAMunsortedUseSoloTmp`. When the two-pass path is selected, a deterministic temp file (`<prefix>Aligned.out.unsorted.solo.tmp`) is opened once and the BGZF handle is left null to prevent accidental writes.
- `source/Parameters.h`: adds `bool outBAMunsortedUseSoloTmp` to propagate the decision across components.
- `source/ParametersSolo.h/.cpp`:
  - Parses `--soloAddTagsToUnsorted` (`addTagsToUnsorted`/`addTagsToUnsortedStr`).
  - Parses `--soloWriteTagTable` (`writeTagTableEnabled`, `writeTagTablePath`); instantiates a shared `BAMTagBuffer` when enabled.
  - Extends validation so CB/UB SAM attributes are allowed when either sorted BAM, two-pass unsorted output, or tag-table export is selected. Misconfigurations emit fatal errors with remediation guidance.
  - Enforces feature requirements (`Gene`, `GeneFull`, `GeneFull_Ex50pAS`, `GeneFull_ExonOverIntron`) when a tag table is requested to guarantee `readInfo` coverage.
  - Destructor cleans up `bamTagBuffer` to avoid leaks in long-lived processes.

## Pass 1 – Temporary Capture (`BAMoutputSoloTmp`)
- Implemented in `source/BAMoutputSoloTmp.{h,cpp}` and wired via `source/ReadAlignChunk.cpp`.
- When `outBAMunsortedUseSoloTmp` is true, each mapping chunk constructs a `BAMoutputSoloTmp` pointing at the shared temp stream; optional access to `BAMTagBuffer` preserves metadata for tag tables.
- Records are stored as `[uint32 payload_len][payload][uint64 trailer]`, where the trailer encodes `iReadAll` and `recordIndex` for downstream lookup.
- Bounds checks (`BAMoutput_oneAlignMaxBytes`) and open-state validation guard against truncated streams. `flush()`/`close()` are exposed for symmetry with other writers.

## Pass 2 – Replay and CB/UB Injection (`BAMunsortedAddSoloTags`)
- Located in `source/BAMunsortedAddSoloTags.{h,cpp}` and invoked from `source/STAR.cpp` after `soloMain.processAndOutput()`.
- Steps per record:
  1. Read header (`size_with_len`) and payload from the temp stream; abort if the size exceeds `BAMoutput_oneAlignMaxBytes` or IO fails.
  2. Append the stored trailer, recover `iReadAll`, and validate that it falls within `SoloFeature::readInfo`.
  3. Call `SoloFeature::addBAMtags`, writing the augmented record to the final BGZF stream (`bgzf_write`).
  4. Monitor size growth against `BAM_ATTR_MaxSize`; log and exit on overflow.
- Debug instrumentation (`STAR_DEBUG_TAG`) decorates buffer allocation, size changes, and readInfo mismatches with timestamps.
- On success the temp file is removed; inability to delete is treated as a warning.

## Tag Table Infrastructure
- `source/BAMTagBuffer.{h,cpp}`: thread-safe accumulator of `BAMTagEntry` (`recordIndex`, `iReadAll`). Internally reserves capacity, guards against `uint32_t` overflow, and exposes `writeTagBinary()` to emit a compact stream.
- `source/BAMTagBinaryWriter.{h,cpp}`: packs `(status_bit, cb_index, umi_packed)` into a header + byte-aligned record stream. Helper utilities (`writeLittleEndian64`, `packAndWrite`, `integerLog2ceil`) guarantee portability.
- `source/SoloFeature_writeTagTable.cpp`: orchestrates export at the end of Solo processing. Determines bit-widths, logs start/finish, hands ownership to `BAMTagBuffer`, and clears the buffer to release memory.
- `source/ReadAlignChunk.cpp`: ensures that tag-table-only runs still route unsorted BAM writes through `BAMoutput` while feeding metadata into `BAMTagBuffer` for later export.

## Safety & Diagnostics
- All new allocations are wrapped with `try/catch` to convert `std::bad_alloc` into actionable error messages.
- Parameter mismatches (unsupported feature sets, missing CB/UB requests) produce fatal errors with explicit "SOLUTION" guidance for parity with upstream messaging.
- `STAR_DEBUG_TAG` toggles high-granularity logging across `BAMTagBuffer` and `BAMunsortedAddSoloTags`, including record counts, monotonicity checks, and final file size validation.

## Testing Summary
- **Deterministic harnesses:**
  - `new/tests/emit_test.sh` – validates binary tag stream and unsorted BAM parity across new flags.
  - `new/tests/integration_test.sh` – end-to-end runs combining CB/UB patching with Solo matrices.
  - `new/tests/tag_test.sh`, `test_synthetic.sh`, `test_trailer_fix.sh` – cover edge cases (multimappers, trailer integrity, synthetic reads).
  - `new/tests/mem_test_tags.sh` – builds with `ASAN=1` and replays stress datasets to confirm absence of leaks/overflows.
- **Artifacts:** inputs reside in `new/testing/data/`; expected outputs and debug captures are in `new/testing/outputs/`.
- **Manual workflows:** documented in `new/scripts/runSTAR*.sh` and `new/docs/two_pass_unsorted_usage.md` for production and debug scenarios.

## Upgrade Considerations
- Users migrating from upstream STAR must update run scripts to toggle `--soloAddTagsToUnsorted` when requesting unsorted BAM output with CB/UB tags.
- Consumers of the binary tag table should decode the 32-byte header (`statusBits`, `cbBits`, `umiBits`, `recordCount`) before iterating records.
- Temp file placement inherits `--outFileNamePrefix`; ensure the filesystem can accommodate the unsorted payload plus trailer overhead when the flag is active.
