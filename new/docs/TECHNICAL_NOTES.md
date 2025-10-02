# Technical Notes – STARsolo Unsorted BAM & Tag Table Extensions

## Scope
This document captures the code changes introduced in the fork relative to upstream STAR 2.7.11b, focusing on unsorted BAM CB/UB injection and tag-table export. File paths refer to the repository root.

## Parameter & Configuration Flow
- `source/Parameters.cpp`: postpones unsorted BAM writer selection until after `pSolo.initialize()`, storing the decision in `outBAMunsortedUseSoloTmp`. When the two-pass path is selected, a deterministic temp file (`<prefix>Aligned.out.unsorted.solo.tmp`) is opened once and the BGZF handle is left null to prevent accidental writes.
- `source/Parameters.h`: adds `bool outBAMunsortedUseSoloTmp` to propagate the decision across components.
- `source/ParametersSolo.h/.cpp`:
  - Parses `--soloAddTagsToUnsorted` (`addTagsToUnsorted`/`addTagsToUnsortedStr`).
  - Parses `--soloWriteTagTable` (`writeTagTableEnabled`, `writeTagTablePath`); instantiates a shared `BAMTagBuffer` when enabled.
  - Parses `--soloSkipProcessing` (`skipProcessingStr`, `skipProcessing`) with validation and summary echo.
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

## Skip Mode – Minimal Post-Processing Path (`--soloSkipProcessing yes`)
### Goal
Run the standard mapping pipeline, then bypass Solo counting/matrix generation while still finalizing per-read artifacts (CB/UB tag injection into the unsorted BAM and optional binary tag table).

### Control Flow
- Entry point: `Solo::processAndOutput()`
  - Early in this function, after barcode statistics, check `pSolo.skipProcessing`.
  - If true, invoke `finalizeMinimalSoloOutputs()` and return.

### `finalizeMinimalSoloOutputs()` Responsibilities
1. Ensure that the feature referenced by `pSolo.samAttrFeature` (typically `GeneFull`) has been constructed.
2. Prepare `soloFeat[feature]->readInfo` so that `BAMunsortedAddSoloTags` can safely call `SoloFeature::addBAMtags` for every record.
3. If `pSolo.writeTagTableEnabled`, call `writeTagTableIfRequested(false)` to flush the tag table without allocating matrices.
4. Emit a concise log line to stdout and `Log.out`: "..... skipping Solo counting (soloSkipProcessing=yes)".

### Algorithmic Adjustments
- `SoloFeature::processRecords()` splits into two paths:
  - Full path (default): performs counting, deduplication, cell filtering, matrix writes, stats.
  - Read-info-only path (skip mode): invokes a lightweight helper (e.g., `prepareReadInfoOnly()`) that does the minimal work needed to fill `readInfo` without allocating matrix structures such as `countCellGeneUMI` and `countMatMult`.
- Allocation guards: expensive vectors and matrices are guarded with `if (!pSolo.skipProcessing)` to avoid unnecessary memory use.
- Unsorted BAM tag injection remains unchanged; the presence of `readInfo` is the only requirement.

### Outputs and Invariants
- `Solo.out/` exists in both modes for compatibility with downstream tools.
- `Solo.out/GeneFull`:
  - Baseline: contains matrices and stats.
  - Skip mode: directory exists but is empty.
- `Aligned.out.bam` and `Aligned.out.cb_ub.bin` are bitwise identical between baseline and skip runs given the same inputs and parameters.

### Validation
- `new/tests/skip_processing_test.sh` executes baseline and skip mode with identical parameters except `--soloSkipProcessing yes` and separate prefixes. It asserts:
  - Presence of `Solo.out` in both outputs, `GeneFull` non-empty for baseline and empty for skip.
  - `diff` parity of `Aligned.out.bam` and `Aligned.out.cb_ub.bin` across runs.
  - Equal BAM line counts (`samtools view | wc -l`).
  - ZG tag presence and non-emptiness in both BAMs.

### Gotchas & Edge Cases
- If neither `--soloAddTagsToUnsorted` nor `--soloWriteTagTable` is enabled, skip mode will produce no per-read artifacts; emit a warning.
- Feature requirements: ZG/ZX population depends on `GeneFull` being present in `--soloFeatures`. Enforce or warn accordingly.
- Incompatible configurations: modes that fundamentally rely on matrices (e.g., certain specialized Solo features) should be rejected while in skip mode with actionable error messages.
- Memory accounting: ensure that matrix-related `clearLarge()` calls are not invoked in skip mode (they should be unreachable) to avoid double-free patterns.
- Determinism: preserving bitwise identity requires stable ordering; audit any logging or non-deterministic iteration that might affect tag-table emission.


## Tag Table Infrastructure
- `source/BAMTagBuffer.{h,cpp}`: thread-safe accumulator of `BAMTagEntry` (`recordIndex`, `iReadAll`). Internally reserves capacity, guards against `uint32_t` overflow, and exposes `writeTagBinary()` to emit a compact stream.
- `source/BAMTagBinaryWriter.{h,cpp}`: packs `(status_bit, cb_index, umi_packed)` into a header + byte-aligned record stream. Helper utilities (`writeLittleEndian64`, `packAndWrite`, `integerLog2ceil`) guarantee portability.
- `source/SoloFeature_writeTagTable.cpp`: orchestrates export at the end of Solo processing. Determines bit-widths, logs start/finish, hands ownership to `BAMTagBuffer`, and clears the buffer to release memory.
- `source/ReadAlignChunk.cpp`: ensures that tag-table-only runs still route unsorted BAM writes through `BAMoutput` while feeding metadata into `BAMTagBuffer` for later export.

### ReadInfo Sink Abstraction
- `source/SoloReadFeature_inputRecords.cpp` now offers an overload that accepts a sink callback to record `(readId, cbIdx, umi, status)` directly. This removes the need for a temporary legacy `readInfo` buffer when `SOLO_USE_PACKED_READINFO` is enabled.
- Call sites in `SoloFeature_countCBgeneUMI.cpp` and `SoloFeature_prepareReadInfoOnly.cpp` pass a lambda that routes to `SoloFeature::recordReadInfo`, maintaining parity and reducing memory overhead.

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

<<<<<<< HEAD
## Custom Gene Annotation Tags (ZG/ZX)
### Overview
- Two new BAM tags that provide enhanced gene annotation beyond standard GX/GN tags:
  - **ZG**: Comma-separated list of Ensembl gene IDs for all overlapping genes
  - **ZX**: Genomic overlap classification (`exonic`, `intronic`, `none`, `spanning`)
- Implementation: `source/ZGZXTags.{h,cpp}` with emission in `source/ReadAlign_alignBAM.cpp`
- Uses `GeneFull` feature annotations to capture both exonic and intronic matches

### Tag Emission Logic
- **ZG Tag Format:**
  - Empty overlap: `ZG:Z:-`
  - Single gene: `ZG:Z:ENSG00000103024`
  - Multiple genes: `ZG:Z:ENSG00000103024,ENSG00000171824`
  - Implementation: `ZGZXTags::formatZGTag()` iterates `annFeat.fSet` and maps indices to `Transcriptome::geID`

- **ZX Tag Mapping:**
  - Derived from `ReadAnnotFeature::ovType` enum values
  - Collapses internal enums to user-facing categories:
    - `exonic`: includes `exonic`, `exonicAS`, `exonic50p`, `exonic50pAS`
    - `intronic`: includes `intronic`, `intronicAS`
    - `none`: includes `none`, `intergenic`
    - `spanning`: fallback for other types
  - Implementation: `ZGZXTags::formatZXTag()` and `overlapTypeToString()`

### Integration Points
- **Modified Files:**
  - `source/IncludeDefine.h`: Added `ATTR_ZG` and `ATTR_ZX` attribute IDs
  - `source/Parameters.h`: Extended `outSAMattrPresent` struct with ZG/ZX fields
  - `source/Parameters_samAttributes.cpp`: Added ZG/ZX parsing logic
  - `source/ReadAlign_alignBAM.cpp`: Emits tags using `GeneFull` annotations
  - `source/ReadAnnotations.h`: Added missing `#include <set>` for compilation
  - `source/Makefile`: Included `ZGZXTags.o` in build

- **Data Flow:**
  1. `Transcriptome::geneFullAlignOverlap()` populates `annFeat.fSet` with gene indices
  2. `Transcriptome_classifyAlign.cpp` sets `annFeat.ovType` based on overlap analysis
  3. `ReadAlign_alignBAM.cpp` calls formatting functions during BAM record construction
  4. Tags are appended to BAM attributes via `bamAttrArrayWrite()`

### Configuration Requirements
- Must include `ZG ZX` in `--outSAMattributes`
- Must include `GeneFull` in `--soloFeatures` (required for gene set population)
- Requires `--outSAMtype BAM` (ZG/ZX are BAM-only, not emitted in SAM format)
- Recommended: `--soloStrand Unstranded` to maximize gene detection

### Validation and Testing
- **Test Script:** `new/scripts/validate_zg_zx.py` verifies tag correctness and format compliance
- **Test Results:** Provides 28% more gene annotations compared to GX tags with 99.8% read coverage
- **Example Output:** Successfully tested on production datasets with genes ENSG00000103024, ENSG00000171824
- Comprehensive documentation: `new/docs/ZG_ZX_Implementation_Summary.md`

### Design Rationale
- **GeneFull vs Gene:** Uses `GeneFull` to capture intronic matches, not just transcript-based exonic hits
- **Multiple Gene Support:** Comma-separated format allows reporting of overlapping genes (e.g., nested/overlapping loci)
- **BAM-Only:** Complex tag format makes SAM output impractical; BAM binary representation is more efficient
- **Strand Handling:** Unstranded mode avoids empty annotations from strand mismatch issues

=======
>>>>>>> da05a276c7ca890005f7d6cfe643a08adb8418ba
## Upgrade Considerations
- Users migrating from upstream STAR must update run scripts to toggle `--soloAddTagsToUnsorted` when requesting unsorted BAM output with CB/UB tags.
- Consumers of the binary tag table should decode the 32-byte header (`statusBits`, `cbBits`, `umiBits`, `recordCount`) before iterating records.
- Temp file placement inherits `--outFileNamePrefix`; ensure the filesystem can accommodate the unsorted payload plus trailer overhead when the flag is active.
<<<<<<< HEAD
- For ZG/ZX tags: add `ZG ZX` to `--outSAMattributes` and ensure `GeneFull` is in `--soloFeatures`.
=======
>>>>>>> da05a276c7ca890005f7d6cfe643a08adb8418ba
