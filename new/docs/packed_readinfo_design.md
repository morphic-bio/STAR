# Design: 8‑Byte Packed Read Info for CB/UMI (STAR Solo Fork)

## Status

Phase: Scaffolded (initial code skeletons added under compile flag `SOLO_USE_PACKED_READINFO`)  
Goal: Replace legacy `readInfoStruct` (16-byte padded) with packed 64-bit word per read.  
Branch: `design/packed-readinfo-8byte`  

---

## 1. Motivation

Current per-read structure:
```c++
struct readInfoStruct {
    uint64 cb;
    uint32 umi;
}; // sizeof likely 16 due to padding
```
Memory grows linearly with reads; at large scale (≥1B reads) this exceeds practical RAM.

Packed approach:
- Store CB index, packed UMI, status in a single 64-bit word.
- Reduce memory by ~50%.
- Normalize data access for tag table and BAM augmentation.
- Provide future extensibility (spare bits, extra status codes).

---

## 2. Bit Layout

Logical (LSB → MSB):
```
| UMI (umiBits) | CB index (cbBits) | Status (3 bits) | Spare (remaining up to bit 63) |
```

Definitions:
- `umiBits = 2 * umiLength` (≤ 32)
- `cbBits = ceil(log2(whitelistSize + 1))` (+1 sentinel)
- `statusBits = 3` (currently)
- `spareBits = 64 - (umiBits + cbBits + 3)`

Example (umiLength=12, wl≈737k):  
`umiBits=24`, `cbBits=20`, `statusBits=3` → used=47, spare=17.

Status codes (initial):
| Value | Meaning |
|-------|---------|
| 0 | Missing / invalid CB |
| 1 | Valid CB + valid UMI |
| 2 | Valid CB + invalid UMI (e.g. homopolymer, N) |
| 3 | Reserved (future) |
| 4–7 | Unused (available) |

---

## 3. Scaffold Reference Map

| File | Purpose | Compile Flag Path |
|------|---------|-------------------|
| `source/PackedReadInfo.h` | Declares `PackedReadInfo` + layout | Always included if `SOLO_USE_PACKED_READINFO` |
| `source/PackedReadInfo.cpp` | Implements packing logic | Same |
| `source/SoloFeature.h` | Adds conditional member `PackedReadInfo packedReadInfo` | Guarded |
| `source/SoloFeature_writeTagTable.cpp` | Switches between legacy and packed binary write | Guarded |
| `source/BAMTagBuffer.h` | Adds `writeTagBinaryPacked` overload | Guarded |
| `source/BAMTagBuffer.cpp` | Implements packed write path | Guarded |
| `source/SoloFeature_addBAMtags.cpp` | Scaffold for retrieving packed CB/UMI and injecting tags | Guarded |
| `source/Makefile` | Must list `PackedReadInfo.o` in `OBJECTS` | Manual change |
| `new/docs/packed_readinfo_design.md` | This design doc | N/A |

---

## 3A. Scaffold Gap Review (2025-02-XX)

- `SoloFeature_countCBgeneUMI.cpp`: still unconditionally resizes `readInfo`; add a compile-time branch that calls `initPackedReadInfo(nReadsInput)` when the packed flag is enabled and avoid allocating the legacy vector.
- `SoloFeature_prepareReadInfoOnly.cpp`: allocates and populates `readInfo`; mirror logic with `packedReadInfo` (including status defaults) so skip-processing mode works with the packed path.
- `SoloFeature_collapseUMIall.cpp`: multiple `readInfo.size()>0` guards need to become compile-time branches; replace direct `readInfo[...]` writes with `packedReadInfo.set()` / `setStatus()` and a helper that encapsulates post-correction UMI updates.
- `SoloFeature_countVelocyto.cpp`: pulls CB/UMI from the gene feature via `->readInfo`; introduce an accessor that returns either legacy or packed data so velocyto counting stays functional under the flag.
- `SoloFeature::finalizeTagTableFromReadInfo()`: currently empty; either delete and inline the remaining work in the writer or implement status reconciliation (e.g. downgrade multi-gene UMIs to status ≠1 before writing).
- `BAMTagBuffer::writeTagBinaryPacked()`: double-check CB sentinel handling (`cbIdx` should mirror the legacy `+1`), and record a TODO for preserving >1 status bit in the file header if needed.
- `PackedReadInfo::maskUMI()`: guard against `umiBits == 0` (current `init` throws) and document behaviour for assays that disable UMIs entirely.
- Audit other modules that touch `readInfo` (e.g. `SoloFilteredCells`, per-feature skip paths) and route them through a shared helper instead of copy-pasted `#ifdef` blocks.

## 4. Initialization Flow

1. Determine total reads: existing pipeline already counts / increments `iReadAll`.
2. After total read capacity is known (counting path) or when `SoloFeature::prepareReadInfoOnly()` is invoked (skip-processing path), branch on the compile flag:
   - Legacy: `readInfo.resize(nReadsInput, {(uint64)-1, (uint32)-1});`
   - Packed: `packedReadInfo.init(nReadsInput, pSolo.cbWLstr.size(), pSolo.umiL);`
3. During barcode + UMI parsing:
   ```c++
   // Compute cbIdx (0 if none), umiPacked (0 if invalid), status
   packedReadInfo.set(iReadAll, cbIdx, umiPacked, status);
   ```
4. Status codes assigned once (update with `setStatus()` only if late rejection occurs).
5. Add `SoloFeature::recordReadInfo(readId, cbIdx, umiPacked, status)` (thin wrapper) so collapse, skip mode, and velocyto paths share the exact packing/transcoding code.

---

## 5. API Summary

```c++
void init(nReads, whitelistSize, umiLength);
void set(readId, cbIdx, umiPacked, status);
uint32_t getCB(readId) const;
uint32_t getUMI(readId) const;
uint8_t  getStatus(readId) const;
void setStatus(readId, status);
```

Helper internal methods (`maskUMI`, `maskCB`, `pack`, `unpack`) are inline for performance.

---

## 6. Integration Points

| Location | Action |
|----------|--------|
| Legacy `readInfo` writes | Replace with `packedReadInfo.set()` |
| BAM tag injection | Use getters; reconstruct strings |
| Tag table binary write | Use `writeTagBinaryPacked()` path |
| UMI collapse (if mutating) | Ensure post-collapse CB/UMI stable before final usage |
| Memory accounting | Replace references to old vector size in logs (optional) |
| Skip-processing fast path | Call shared `recordReadInfo()` helper so `prepareReadInfoOnly()` populates packed storage |
| Velocyto counting | Read CB/UMI/status via gene feature accessor that dispatches to `packedReadInfo` |
| Downstream analytics | Any module pulling from `soloFeatAll[gene]->readInfo` must migrate to helper API |

---

## 7. BAM Tag Injection (Scaffold Notes)

Current scaffold (`SoloFeature_addBAMtags.cpp`) includes a placeholder UMI decode:

```c++
for (int i = pSolo.umiL - 1; i >= 0; --i) {
    uint8_t b = tmp & 0x3;
    ...
}
std::reverse(umiStr.begin(), umiStr.end());
```

If your existing packing order differs (e.g., first base in high bits), adjust to maintain parity with legacy path.

Status handling: only emit CB/UB strings when `status == 1`; for other statuses, emit `"-"` (legacy behaviour) and consider dropping UB entirely if UMI packer stored zero.
Ensure CB index parity: legacy code stores WL index; packed path must mirror that and guard against OOB indexes before dereferencing `pSolo.cbWLstr`.

Inject tags with your existing helper (example placeholders):

```c++
addBAMtagString(bamRec, bamSize, "CB", cbStr.c_str());
addBAMtagString(bamRec, bamSize, "UB", umiStr.c_str());
```

---

## 8. Binary Tag Table Writer

Legacy path uses:
```c++
writeTagBinary(path, readInfo, cbBits, umiBits);
```

Packed path:
```c++
writeTagBinaryPacked(path, packedReadInfo, cbBits, umiBits);
```

Status mapping to 1-bit:
```c++
statusBit = (packedStatus == 1 ? 1 : 0);
```
(Extend header later if multi-level statuses need to be preserved.)

Sentinel alignment:
- Legacy writer adds `+1` to CB index before serialising (0 = sentinel). Ensure packed path performs the same offset or document the change and bump the file format version.
- Add TODO: if future statuses require >1 bit, version the `TagStreamHeader` and persist the full status when header advertises it.

---

## 9. Compile Flag Strategy

| Flag | Behavior |
|------|----------|
| (unset) | Legacy `readInfoStruct` used |
| `-DSOLO_USE_PACKED_READINFO` | Packed implementation compiled; legacy retained in conditional blocks |

Migration Phases:
1. Dual-mode (default: legacy)
2. Default flipped (packed by default)
3. Legacy removed

---

## 10. Test Matrix

| Scenario | Goal | Execution | Notes |
|----------|------|-----------|-------|
| BAM parity (10k reads) | Detect functional regressions | Run pipeline twice (`make clean && make`, toggle `-DSOLO_USE_PACKED_READINFO`), compare CB/UB tags | Expect byte-for-byte equality in BAM tags; mismatch ⇒ inspect decoding order |
| Binary tag stream diff | Validate sidecar writer | Compare `*.solo/CB*.bin` between legacy and packed builds | Header may differ only in metadata version; record bytes must match after CB `+1` offset applied |
| Invalid UMI bases | Ensure status codes propagate | Inject reads with `N` in UMI, confirm `status==2`, UMI bits zero | Update test to assert `UB` tag becomes `-` |
| Missing CB path | Confirm sentinel behaviour | Use reads outside whitelist; expect `status==0`, CB/UMI zero | Verify packed writer emits 0 sentinel and does not segfault on WL lookup |
| Whitelist edge (2^k boundary) | Exercise cbBits growth | Run with synthetic WL one below / above power of two | Ensure `log2ceil` matches legacy; add unit test for `PackedReadInfo::init` |
| UMI length extremes (4/12/16/18) | Check bit budget | Unit test `PackedReadInfo::init` and `pack/unpack` for supported lengths | For >16 (⇒ 32 bits) confirm guard throws with clear message |
| Skip-processing pipeline | Verify read info population w/out counting | Execute `runSTAR.sh` with `--soloSkipProcessing yes` under both builds | Ensure `writeTagTableIfRequested` works with packed storage |
| Velocyto counting | Ensure downstream modules read packed data | Run velocyto-enabled config; compare counts | Requires accessor bridging to packed data |
| Memory footprint | Confirm RAM reduction | `/usr/bin/time -v` around STAR run, capture `Max resident set size` | Expect ~50% reduction for read info component |
| Bit audit (fuzz) | Catch packing bugs | Add googletest/standalone unit: random `cbIdx/umi/status`, `pack→unpack` stable | Include negative tests for status overflow |

Suggested automated validation snippet:
```bash
# After running both builds from a clean tree
python compare_cb_ub.py legacy.bam packed.bam  # script parses CB/UB tags, diffs sets
python compare_cb_ub.py legacy_tag.bin packed_tag.bin --binary
```

---

## 10A. Test Execution Checklist

1. Build legacy baseline (`make clean && make STAR`) and archive `Aligned.out.bam` plus tag binaries.
2. Rebuild with `CXXFLAGS+=' -DSOLO_USE_PACKED_READINFO'` (or export `STAR_USE_PACKED=1` wrapper) and rerun identical input.
3. Run comparison scripts:
   - `python scripts/verify_packed_vs_legacy.py legacy.bam packed.bam`
   - `python scripts/compare_tag_bin.py legacy_tag.bin packed_tag.bin`
4. Execute targeted unit tests (once added) via `ctest -R PackedReadInfo` or standalone binary.
5. Capture `/usr/bin/time -v` output for both runs; stash in `new/docs/perf/packed_readinfo_rss.txt` for trend tracking.
6. Record results in `DEBUG_RUN_STATUS.md` with configuration hashes (commit + compile flag) for reproducibility.

---

## 11. Failure Modes & Guards

| Failure | Guard |
|---------|-------|
| `usedBits > 64` | Throw in `init()` |
| `umiBits > 32` (UMI too long) | Throw (or future fallback) |
| `cbBits > 32` (WL huge) | Throw; instruct fallback |
| `readId >= data.size()` | Debug assert (optionally) |
| Status code overflow | Mask with `(1<<statusBits)-1` |
| CB sentinel mismatch | Unit test ensures `getCB` returns WL index; writer applies `+1` |
| Skip-processing bypasses init | Guard `prepareReadInfoOnly()` to call helper under compile flag |

---

## 12. Performance Considerations

- All operations O(1) bit shifts/masks.
- Inline decoding negligible relative to alignment + parsing cost.
- Sort cost unchanged (tag buffer still sorted by `recordIndex`).
- Reduced memory footprint improves cache residency.

---

## 13. Extension Slots (Spare Bits)

Immediate reserved for:
- Future: UMI collision collapse state
- Gene multiplicity flag
- Short hash for rapid duplicate detection

Document high-bit usage before assignment if extended.

---

## 14. Decommissioning Legacy Path

Checklist when ready:
- Remove `#ifndef SOLO_USE_PACKED_READINFO` branches.
- Delete `readInfoStruct` (if nowhere else referenced).
- Update `TECHNICAL_NOTES.md`.
- Remove compile flag from default build pipeline.
- Add regression test baseline hash updated.

---

## 15. Diff Example (Call Site Migration)

```diff
- readInfo[iRead].cb  = cbIndex;
- readInfo[iRead].umi = umiPacked;
+ packedReadInfo.set(iRead, cbIndex, umiPacked, statusCode);
```

```diff
- const readInfoStruct &R = readInfo[entry.readId];
- uint64_t cb = R.cb;
- uint32_t umi = R.umi;
+ uint64_t cb = packedReadInfo.getCB(entry.readId);
+ uint32_t umi = packedReadInfo.getUMI(entry.readId);
+ uint8_t status = packedReadInfo.getStatus(entry.readId);
```

---

## 16. Implementation Refinement Checklist

| Step | Action | Done |
|------|--------|------|
| 1 | Add new files (PackedReadInfo.h/.cpp) | ☑ |
| 2 | Add object to Makefile | ☑ |
| 3 | Guard legacy vs packed in `SoloFeature.h` | ☑ |
| 4 | Initialize packed in constructor / setup | ☑ |
| 5 | Implement `recordReadInfo()` helper shared by legacy/packed | ☑ |
| 6 | Wire skip-processing (`prepareReadInfoOnly`) through helper | ☑ |
| 7 | Replace collapse paths with `.set()` / `setStatus()` | ☑ |
| 8 | Provide velocyto accessor for packed data | ☑ |
| 9 | Adjust BAM tag injection | ☑ |
| 10 | Adjust tag binary writer | ☑ |
| 11 | Run parity tests | ☐ |
| 12 | Measure RSS difference | ☐ |
| 13 | Add automated test (CI) | ☐ |
| 14 | Flip default (enable packed w/out flag) | ☐ |
| 15 | Remove legacy code | ☐ |

---

## 17. Rollback Plan

If issues arise with packed path:
1. Rebuild without `-DSOLO_USE_PACKED_READINFO`.
2. Compare outputs to isolate divergence (likely UMI decode ordering or status mapping).
3. Re-enable packed after fix.

---

## 18. Documentation Addition (for TECHNICAL_NOTES.md)

Add snippet:

```
READ INFO PACKING:
Each read's corrected CB index, packed UMI, and status are stored in a 64-bit word:
  [ UMI (2*L bits) | CB (ceil(log2(WL+1))) | Status (3) | Spare ]
Legacy readInfoStruct removed after version X.Y.Z.
```

---

## 19. Future (Optional) Enhancements

| Idea | Description |
|------|-------------|
| 6-byte variant | Further compress when `umiBits+cbBits+statusBits ≤ 48` |
| Extended status | Use spare bits for multi-class failure labels |
| On-the-fly direct sidecar | Emit alignment records as they come (skip post-sort write) |
| Memory-mapped tag file | Avoid large resident buffers during giant runs |

---

## 20. Conclusion

The scaffolds are integrated behind a compile flag and referenced here. Implementation should proceed by migrating call sites with the checklist, validating parity, and then removing the legacy structure.

---

## 21. Appendix: Quick Verification Script (Example)

```python
# scripts/verify_packed_vs_legacy.py
import pysam, sys

def collect(fn):
    d={}
    with pysam.AlignmentFile(fn, 'rb') as bam:
        for i,rec in enumerate(bam):
            cb = rec.get_tag('CB') if rec.has_tag('CB') else None
            ub = rec.get_tag('UB') if rec.has_tag('UB') else None
            d[i]=(cb,ub)
    return d

legacy = collect(sys.argv[1])
packed = collect(sys.argv[2])

missing=0
mismatch=0
for k,v in legacy.items():
    if k not in packed: missing+=1
    elif packed[k]!=v: mismatch+=1

print("Missing:", missing, "Mismatch:", mismatch)
assert mismatch==0
```

Run:
```
python scripts/verify_packed_vs_legacy.py legacy.bam packed.bam
```

---

*End of Document*
