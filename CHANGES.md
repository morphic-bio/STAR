# Changes Since the STAR 2.7.11b Fork

This fork extends upstream STAR 2.7.11b with tighter support for STARsolo runs that emit unsorted BAM files and/or standalone barcode tables.

## Feature Highlights

### 1. Two-Pass CB/UB Injection for Unsorted BAM (`--soloAddTagsToUnsorted`)
- Introduces `--soloAddTagsToUnsorted` (`yes`|`no`, default `no`). When enabled alongside `--outSAMtype BAM Unsorted`, STAR captures unsorted alignments in a temporary stream during pass 1 and rewrites the final BAM after Solo correction.
- Pass 1 buffering is handled by `BAMoutputSoloTmp` (`source/BAMoutputSoloTmp.{h,cpp}`), wired through `ReadAlignChunk` to collect `[size][payload][iReadAll]` triplets per alignment.
- Pass 2 replay is performed by `BAMunsortedAddSoloTags` (`source/BAMunsortedAddSoloTags.{h,cpp}`), which rehydrates the stream, calls `SoloFeature::addBAMtags`, validates record sizes, writes the final BAM, and removes the temporary file.
- Parameters guard against invalid configurations (e.g., CB/UB attributes with unsorted BAM unless either the new flag or tag-table export is enabled) and provide explicit error messaging.

**Usage Example:**
```bash
STAR \
  --genomeDir /path/to/genome \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloFeatures Gene GeneFull \
  --outSAMtype BAM Unsorted \
  --soloAddTagsToUnsorted yes \
  --outFileNamePrefix /mnt/run1/
```

### 2. Binary Tag Table Export (`--soloWriteTagTable`)
- Adds `--soloWriteTagTable` with accepted values `None` (default), `Default`, or a custom path. When enabled, corrected CB/UB assignments are streamed into a shared `BAMTagBuffer` and emitted as a compact binary sidecar (`Aligned.out.cb_ub.bin` by default).
- `BAMTagBinaryWriter` packs status, CB index, and packed UMI into a header+record stream; `SoloFeature::writeTagTableIfRequested` orchestrates emission, logging, and buffer lifecycle management.
- The export pathway can be used with or without unsorted BAM patching, allowing lightweight audits or secondary processing pipelines.

**Usage Example (sidecar only):**
```bash
STAR \
  --genomeDir /path/to/genome \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloFeatures Gene \
  --outSAMtype BAM Unsorted \
  --soloWriteTagTable Default \
  --soloAddTagsToUnsorted no
```

**Usage Example (sidecar + patched BAM):**
```bash
STAR \
  --genomeDir /path/to/genome \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --soloType CB_UMI_Simple \
  --outSAMtype BAM Unsorted \
  --soloAddTagsToUnsorted yes \
  --soloWriteTagTable /tmp/run.cb_ub.bin
```

<<<<<<< HEAD
### 3. Custom Gene Annotation Tags (ZG/ZX)
- Introduces two new BAM tags that enhance gene annotation capabilities beyond standard GX/GN tags:
  - **ZG**: Gene set tag containing comma-separated Ensembl gene IDs for all genes overlapping the read
  - **ZX**: Overlap classification tag indicating genomic region type (`exonic`, `intronic`, `none`, `spanning`)
- Uses `GeneFull` feature annotations to capture both exonic and intronic matches
- Provides 28% more gene annotations compared to standard GX tags with 99.8% read coverage
- Implemented in `ZGZXTags.{h,cpp}` with proper SAM tag formatting

**ZG Tag Format:**
- Single gene: `ZG:Z:ENSG00000103024`
- Multiple genes: `ZG:Z:ENSG00000103024,ENSG00000171824`
- No genes: `ZG:Z:-`

**ZX Tag Values:**
- `exonic`: Read overlaps exons (includes antisense and ≥50% exonic cases)
- `intronic`: Read overlaps introns (includes antisense intronic)
- `none`: No gene overlap (intergenic)
- `spanning`: Other/unknown overlap types

**Usage Example:**
```bash
STAR \
  --genomeDir /path/to/genome \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --outSAMattributes NH HI AS nM NM CR CY UR UY GX GN ZG ZX \
  --soloType CB_UMI_Simple \
  --soloFeatures Gene GeneFull \
  --soloStrand Unstranded \
  --outFileNamePrefix output/
```

**Key Requirements:**
- Must include `ZG ZX` in `--outSAMattributes`
- Must include `GeneFull` in `--soloFeatures` 
- Requires `--outSAMtype BAM` (ZG/ZX are BAM-only tags)
- Recommended: `--soloStrand Unstranded` to avoid strand specificity issues

For detailed implementation information, see `new/docs/ZG_ZX_Implementation_Summary.md`.

=======
>>>>>>> da05a276c7ca890005f7d6cfe643a08adb8418ba
### 3. Solo Skip Processing (`--soloSkipProcessing {yes|no}`)
- Adds `--soloSkipProcessing` (default `no`) to bypass Solo counting, UMI deduplication, matrix generation, and cell filtering while preserving per-read artifacts. The mapping pipeline runs to completion, ZG/ZX and CB/UB tags are finalized, and the optional binary tag table is written if enabled.
- Baseline vs Skip expectations:
  - Identical `Aligned.out.bam` and `Aligned.out.cb_ub.bin`
  - `Solo.out/` directory exists in both runs; `Solo.out/GeneFull` contains files in baseline and is empty in skip mode
  - ZG tag fields are populated (non-empty and not `-`) in BAM
- Parameter surfaces:
  - `source/ParametersSolo.h/.cpp`: `skipProcessingStr` and parsed `skipProcessing` boolean with validation and logging
  - `source/parametersDefault`: added `soloSkipProcessing           no`
  - `source/Parameters.cpp`: CLI registration and summary echo
- Control-flow hook:
  - `Solo::processAndOutput()` branches early when `skipProcessing` is true to execute a minimal post-processing path that prepares `readInfo` for tag injection and optional tag-table emission, then returns without matrix allocations
- Safety/validation:
  - Warns when skip mode is enabled without any per-read outputs requested (e.g., neither `--soloAddTagsToUnsorted` nor `--soloWriteTagTable`)
  - Incompatible feature sets (that fundamentally require matrices) are rejected with actionable errors
- Tests:
  - `new/tests/skip_processing_test.sh` runs baseline and skip side-by-side, confirms Solo folder expectations, diffs BAM and BIN parity, checks equal BAM line counts, and validates non-empty ZG fields in both outputs

## Infrastructure Updates
- `Parameters.cpp` defers unsorted BAM opening until Solo initialization has resolved whether the two-pass workflow is required, storing the decision in `outBAMunsortedUseSoloTmp`.
- `ParametersSolo` parses and validates both new flags, manages the shared `BAMTagBuffer`, and extends error handling for CB/UB attribute combinations and feature requirements.
- Debug hooks (gated by `STAR_DEBUG_TAG`) add timestamped diagnostics across the new modules without affecting production builds.
- Safety checks cover buffer allocation failures, oversized records, and invalid read indices during replay.

## Testing & Validation
- Deterministic harnesses (`new/tests/emit_test.sh`, `integration_test.sh`, `mem_test_tags.sh`, `tag_test.sh`, `test_synthetic.sh`, `test_trailer_fix.sh`) were updated to exercise the new flags.
- Golden outputs and large-scale validation artifacts are archived under `new/testing/outputs/` with input data in `new/testing/data/`.
- AddressSanitizer builds (`ASAN=1 make STAR`) and debug-mode runs (`STAR_DEBUG_TAG=1`) verified memory safety and buffer accounting for the new pipelines.

For implementation deep dives, design rationale, and verification matrices, see `new/docs/TECHNICAL_NOTES.md`.
