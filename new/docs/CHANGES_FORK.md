# STARsolo Fork – Changes Since Upstream

This fork extends STARsolo with enhanced CB/UB tag handling for unsorted BAM output. Key changes:

## Added Features

### Two-Pass Unsorted CB/UB Tag Injection
- **New CLI flag – `--soloAddTagsToUnsorted`** (default `no`): when set to `yes` alongside `--outSAMtype BAM Unsorted`, STAR captures alignments in a temporary container, performs Solo correction, then injects the corrected CB/UB tags into the final unsorted BAM.
- **Pass-1 capture**: unsorted BAM output is redirected to a new `BAMoutputSoloTmp` writer that stores `[size][payload][trailer]` triplets containing the alignment record and `iReadAll` identifier.
- **Pass-2 replay**: a new module, `BAMunsortedAddSoloTags`, replays the temporary file, invokes `SoloFeature::addBAMtags`, and writes the final unsorted BAM while preserving tag semantics identical to the coordinate-sorted pipeline.

### Tag Table Export
- **New CLI flag – `--soloWriteTagTable`** (default `None`): exports corrected CB/UB assignments to a sidecar table instead of or in addition to patching them into the BAM.
  - `None` (default): disables tag table export
  - `Default`: writes to `<outFileNamePrefix>Aligned.out.cb_ub.tsv`
  - Custom path: writes to specified file location
- **Sidecar format**: tab-separated table with columns for `bam_record_index`, `iReadAll`, `mate`, `align_idx`, `qname`, `CB`, `UB`, and `status`, strictly aligned with BAM record order.
- **Memory efficient**: bypasses the temporary BAM file entirely when only tag table export is needed.

### Implementation Updates
- **Shared plumbing updates**: parameter parsing, SAM header emission, and STAR's orchestration flow were updated to open/close the appropriate streams based on the new flags. Buffer guardrails were added to prevent oversized records from overflowing fixed scratch buffers.
- **Testing harness**: the `testing/` directory now contains scripts that reproduce the sorted baseline, run the new modes, and compare BAM tags, alignment fields, Solo matrices, and logs to guarantee behavioural parity.

Refer to `docs/two_pass_unsorted_usage.md` for usage instructions and to `docs/two_pass_unsorted_tech.md` for implementation details.
