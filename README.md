# STAR Fork Overview

<<<<<<< HEAD
This repository tracks a focused fork of STAR 2.7.11b that extends STARsolo for unsorted BAM workflows, tag table export, and enhanced gene annotation.
=======
This repository tracks a focused fork of STAR 2.7.11b that extends STARsolo for unsorted BAM workflows and tag table export.
>>>>>>> da05a276c7ca890005f7d6cfe643a08adb8418ba

## Quick Summary
- Two-pass unsorted BAM patching via `--soloAddTagsToUnsorted yes` preserves corrected CB/UB tags when `--outSAMtype BAM Unsorted` is requested.
- Binary tag table export via `--soloWriteTagTable` captures barcode assignments without rewriting the BAM stream.
<<<<<<< HEAD
- Custom gene annotation tags (`ZG`/`ZX`) provide comprehensive gene ID and overlap classification for each aligned read.
=======
>>>>>>> da05a276c7ca890005f7d6cfe643a08adb8418ba
- Shared safety work adds guardrails (buffer bounds, parameter validation, debug hooks) to the new pipelines.

See `CHANGES.md` for the full list of changes, CLI examples, and compatibility notes. Implementation details and test matrices live in `new/docs/TECHNICAL_NOTES.md`.

## Directory Layout
- `source/`, `bin/`, `doc/`, `examples/`, `extras/`, `tools/` – upstream STAR sources and assets (required for compilation).
- `new/` – fork-only material grouped by purpose (`docs`, `plans`, `scripts`, `tests`, `testing`, `src`).

For legacy STAR release history and ancillary documentation, refer to `new/docs/`.
