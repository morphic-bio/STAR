# STAR Fork Overview

This repository tracks a focused fork of STAR 2.7.11b that extends STARsolo for unsorted BAM workflows and tag table export.

## Quick Summary
- Two-pass unsorted BAM patching via `--soloAddTagsToUnsorted yes` preserves corrected CB/UB tags when `--outSAMtype BAM Unsorted` is requested.
- Binary tag table export via `--soloWriteTagTable` captures barcode assignments without rewriting the BAM stream.
- Shared safety work adds guardrails (buffer bounds, parameter validation, debug hooks) to the new pipelines.

See `CHANGES.md` for the full list of changes, CLI examples, and compatibility notes. Implementation details and test matrices live in `new/docs/TECHNICAL_NOTES.md`.

## Directory Layout
- `source/`, `bin/`, `doc/`, `examples/`, `extras/`, `tools/` – upstream STAR sources and assets (required for compilation).
- `new/` – fork-only material grouped by purpose (`docs`, `plans`, `scripts`, `tests`, `testing`, `src`).

For legacy STAR release history and ancillary documentation, refer to `new/docs/`.
