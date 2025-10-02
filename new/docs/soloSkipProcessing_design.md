# Design: `--soloSkipProcessing` (Skip Solo Post-Mapping Counting Phase)

## Summary
Introduce a new CLI flag `--soloSkipProcessing` (default: `no`) to allow STAR (fork) users to:
- Perform alignment, per-read gene overlap classification (`GeneFull`), and emit ZG / ZX tags.
- Optionally inject CB/UB tags into the BAM (`--soloAddTagsToUnsorted yes`) and/or export the CB/UB binary tag table (`--soloWriteTagTable`).
- Skip the memory- and time-intensive Solo counting, feature matrix construction, UMI deduplication, and cell filtering logic executed in `Solo::processAndOutput()`.

This reduces peak RAM after alignment for workflows that require only per-read annotations (ZG/ZX, CB/UB) and not cell-by-feature count matrices.

## Motivation
Current pipeline:
1. Alignment phase builds necessary per-read annotations (enough for ZG/ZX formatting).
2. After alignment, Solo counting allocates and processes large matrices (barcodes × features), causing a secondary memory spike.
3. Many lightweight QC or downstream tools only need:
   - ZG/ZX tags (multi-gene context + overlap classification),
   - CB/UB tags inside BAM or a sidecar binary (audit / reproducibility).
   They do NOT need raw or filtered count matrices.

Providing a skip flag saves wall-clock time and memory, especially on large datasets with high cell or feature cardinality.

## Scope
In scope:
- New flag parsing and validation.
- Early exit path that suppresses Solo feature processing.
- Ensuring ZG/ZX emission continues to work.
- Ensuring CB/UB tag injection and binary tag table finalization still function.

Out of scope:
- Changing existing default behavior.
- Modifying format of existing outputs (ZG/ZX or binary tag table).
- Refactoring upstream (non-fork) code paths.

## Flag Specification
| Flag | Values | Default | Effect |
|------|--------|---------|--------|
| `--soloSkipProcessing` | `yes` / `no` | `no` | When `yes`, bypasses Solo counting, matrix construction, UMI dedup, and cell filtering. |

Behavioral matrix:

| Need ZG/ZX | Need CB/UB in BAM | Need binary tag table | Need matrices | Use `--soloSkipProcessing`? |
|------------|------------------|-----------------------|--------------|-----------------------------|
| Yes | Optional | Optional | No | Yes |
| Yes | Yes | Yes | No | Yes |
| Yes | Yes | Yes | Yes | No |
| No | Yes | Yes | No | Yes (if ZG/ZX not required) |
| Yes | No | No | Yes | No |

## Detailed Behavior
When `--soloSkipProcessing yes`:
1. Alignment proceeds normally.
2. ZG/ZX tags are emitted per alignment if requested in `--outSAMattributes` and `GeneFull` is present in `--soloFeatures`.
3. CB/UB injection (unsorted two-pass) still works if `--soloAddTagsToUnsorted yes`.
4. Binary tag table still written if `--soloWriteTagTable` != `None`.
5. `Solo::processAndOutput()` performs (optionally) minimal barcode stats aggregation and returns early.
6. No matrices, no cell filtering, no UMI dedup passes.
7. Logging notes skip mode clearly.

## Implementation Plan

### 1. Parameters & Data Structures
Add to `ParametersSolo.h`:
```cpp
std::string skipProcessingStr = "no";
bool skipProcessing = false;
```

Register flag in `Parameters.cpp` near other solo parameters:
```cpp
parArray.push_back(new ParameterInfoScalar<string>(-1, -1,
    "soloSkipProcessing", &pSolo.skipProcessingStr));
```

Parse in `ParametersSolo::initialize()`:
```cpp
if (skipProcessingStr == "yes") {
    skipProcessing = true;
} else if (skipProcessingStr == "no") {
    skipProcessing = false;
} else {
    // fatal error with guidance
}
```

### 2. Constructor Adjustment (`Solo::Solo(ReadAlignChunk** ...)
Current logic allocates `SoloFeature` objects unless type is `CB_samTagOut`.
Modify:
```cpp
if (pSolo.type == 0) return;

readBarSum = new SoloReadBarcode(P); // still needed for stats & WL corrections

if (pSolo.type == pSolo.SoloTypes::CB_samTagOut) return;
if (pSolo.skipProcessing) return; // do NOT allocate soloFeat
```

### 3. Early Exit in `Solo::processAndOutput()`
At the top:
```cpp
if (pSolo.type == 0) return;

if (pSolo.skipProcessing) {
    // Optional: gather barcode stats only
    if (pSolo.cbWLyes && readBarSum != nullptr) {
        for (int ii=0; ii<P.runThreadN; ii++) {
            readBarSum->addCounts(*RAchunk[ii]->RA->soloRead->readBar);
            readBarSum->addStats(*RAchunk[ii]->RA->soloRead->readBar);
            delete RAchunk[ii]->RA->soloRead->readBar;
        }
        ofstream *statsStream = &ofstrOpen(
            P.outFileNamePrefix + pSolo.outFileNames[0] + "Barcodes.stats",
            ERROR_OUT, P);
        readBarSum->statsOut(*statsStream);
        statsStream->close();
    }
    if (pSolo.writeTagTableEnabled && pSolo.bamTagBuffer) {
        pSolo.bamTagBuffer->finalize(pSolo.writeTagTablePath);
    }
    *P.inOut->logStdOut << timeMonthDayTime()
        << " ..... skipping Solo counting (soloSkipProcessing=yes)\n";
    *P.inOut->logMain << timeMonthDayTime()
        << " ..... skipping Solo counting (soloSkipProcessing=yes)\n";
    return;
}
```

### 4. ZG/ZX Independence
Ensure gene full overlap calculation runs if:
```cpp
needGeneFull = outSAMattrPresent.ZG || outSAMattrPresent.ZX || GeneFull requested explicitly
```
If current code ties `GeneFull` activation to `soloFeatures`, confirm user still supplies `--soloFeatures GeneFull`; otherwise consider automatically enabling internal geneFull logic when ZG/ZX attributes requested and warn if user omitted it.

### 5. Interaction Validation
In `ParametersSolo::initialize()` after other validations:
- Warn if `skipProcessing && !writeTagTableEnabled && !addTagsToUnsorted`:
  > "WARNING: --soloSkipProcessing=yes but neither CB/UB tag injection nor tag table export is enabled; skipping Solo does not produce matrices or tag sidecar."
- Warn if user also specified `--soloCellFilter*` or cell filtering modes:
  > "NOTE: cell filtering is disabled by --soloSkipProcessing; ignoring related parameters."

### 6. Logging
Add a line to final summary (Log.out):
```
soloSkipProcessing             yes|no
```
Add to help / documentation (`CHANGES.md` and `new/docs/TECHNICAL_NOTES.md`).

### 7. Testing Strategy
| Test | Flags | Expected |
|------|-------|----------|
| Baseline | (no skip) | Matrices present |
| Skip minimal | `--soloSkipProcessing yes --soloFeatures GeneFull --outSAMattributes ZG ZX` | ZG/ZX present; no matrices; lower memory |
| Skip + Tag Table | `--soloSkipProcessing yes --soloWriteTagTable Default` | Binary file exists; no matrices |
| Skip + Tag Injection | `--soloSkipProcessing yes --soloAddTagsToUnsorted yes` | CB/UB tags in BAM; no matrices |
| Conflict Warn | `--soloSkipProcessing yes --quantMode GeneCounts` | GeneCounts unaffected (bulk) or warning if disallowed |
| Error Input | `--soloSkipProcessing maybe` | Fatal parameter error |

Add integration script `new/tests/skip_processing_test.sh` to:
1. Run with/without skip.
2. Parse peak RSS (e.g., `/usr/bin/time -v`).
3. Assert absence of `Solo.out/*/raw/*` matrices in skip mode.
4. Confirm presence of ZG tags: `samtools view | grep -m1 "ZG:Z:"`.

### 8. Performance Expectations
- RAM: Avoids allocating `soloFeat` arrays and large count matrices—savings scale with (#barcodes × #features).
- CPU: Skips counting loops, reducing post-alignment runtime proportionally (often 10–40% for large datasets).
- I/O: Fewer matrix files written.

### 9. Backward Compatibility
Default = `no`; existing workflows unchanged.
Flag gated; only additive risk is ensuring early returns don’t bypass required finalizations (tag table finalization explicitly handled).

### 10. Risks & Mitigations
| Risk | Mitigation |
|------|------------|
| Tag table not finalized | Explicit finalize call in skip path |
| ZG/ZX missing due to mis-specified `--soloFeatures` | Emit warning if ZG/ZX requested but GeneFull absent |
| Users assume matrices exist | Log and WARN clearly when skipping |
| Future code expects `soloFeat` non-null | Guard usages or set `soloFeat=nullptr` and null-check |
| Incomplete WL stats | Optionally gather minimal stats; document difference |

### 11. Rollback Plan
- Removing the flag only requires deleting parsing block, constructor early return, and skip branch in `processAndOutput()`. 
- No persistent format changes, so rollback is trivial.

### 12. Documentation Addendum (CHANGES.md snippet)
```
### Added: --soloSkipProcessing
Skips Solo post-mapping counting and cell filtering while still permitting:
- ZG/ZX BAM tags (requires --soloFeatures GeneFull)
- CB/UB BAM tag injection (--soloAddTagsToUnsorted yes)
- CB/UB tag table export (--soloWriteTagTable)
Reduces peak memory and runtime for per-read annotation workflows.
```

### 13. Example Usage

Minimal ZG/ZX only:
```bash
STAR \
  --genomeDir /ref/star_index \
  --readFilesIn R2.fq.gz R1.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --outSAMattributes NH HI AS nM NM ZG ZX \
  --soloType CB_UMI_Simple \
  --soloFeatures GeneFull \
  --soloSkipProcessing yes \
  --outFileNamePrefix skipDemo/
```

With tag table & CB/UB injection:
```bash
STAR \
  --genomeDir /ref/star_index \
  --readFilesIn R2.fq.gz R1.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --outSAMattributes NH HI AS nM NM CR CY UR UY CB UB ZG ZX \
  --soloType CB_UMI_Simple \
  --soloFeatures GeneFull \
  --soloCBwhitelist whitelist.txt \
  --soloAddTagsToUnsorted yes \
  --soloWriteTagTable Default \
  --soloSkipProcessing yes \
  --outFileNamePrefix skipFull/
```

## Appendix: Quick Diff Summary (Conceptual)
- Parameters.cpp: + registration
- ParametersSolo.h: + members
- ParametersSolo.cpp: + parse/validate, + help text (optional)
- Solo.cpp: + early returns & skip branch
- (Optional) ReadAlign_alignBAM.cpp: ensure ZG/ZX logic is unconditional given attributes

---
## Next Steps
1. Implement code changes.
2. Add tests and doc updates.
3. Benchmark memory before/after on representative dataset.
4. Merge once validated.
