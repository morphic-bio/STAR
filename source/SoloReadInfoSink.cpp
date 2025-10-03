#include "SoloReadInfoSink.h"
#include "SoloFeature.h"
#include "SoloReadInfoLoader.h"

void MinimalSink::onRecord(SoloFeature &feature, const ReadInfoRecord &rec) {
    // Match CountingSink filter: only persist records with valid status and valid feature
    // This keeps both modes in sync and preserves legacy behavior for no-feature reads
    if (rec.status != 1 || rec.featureId == (uint32_t)-1) {
        // Explicitly write sentinel for rejected/no-feature reads to ensure binary writer sees status=0
        feature.recordReadInfo((uint32_t)rec.readId, 0, 0, 0);
        return;
    }
    
    // Pass original values; recordReadInfo will translate per-backend
    feature.recordReadInfo((uint32_t)rec.readId, rec.cbIdx, rec.umi, rec.status);
}

void MinimalSink::finalize(SoloFeature &feature) {
    (void)feature;
}

void CountingSink::onRecord(SoloFeature &feature, const ReadInfoRecord &rec) {
    // Buffer counted records per WL CB for later materialization into rGeneUMI.
    // Only accept good CB/UMI and valid featureId.
    if (rec.status != 1) return;
    if (rec.featureId == (uint32_t)-1) return;
    // Guard: if no whitelist, skip safely
    if (feature.pSolo.cbWLsize==0) return;
    if (perWL.empty()) perWL.assign(feature.pSolo.cbWLsize, {});
    if (rec.cbIdx >= feature.pSolo.cbWLsize) return;
    perWL[rec.cbIdx].push_back(rec);
}

void CountingSink::finalize(SoloFeature &feature) {
    // Materialize buffered triplets into rGeneUMI/rCBp and run collapseUMIall().
    if (perWL.empty()) return;
    // Reset state for this finalize run
    feature.indCB.clear();
    feature.indCBwl.assign(feature.pSolo.cbWLsize, (uint32) -1);
    feature.nCB = 0;
    feature.setRGUStride(feature.pSolo.readIndexYes[feature.featureType] ? 3u : 2u);

    // Build detected CB list in WL order for determinism
    std::vector<uint32_t> cbWLtoICB(perWL.size(), (uint32_t)-1);
    for (uint32_t wl=0; wl<perWL.size(); ++wl) {
        if (!perWL[wl].empty()) {
            cbWLtoICB[wl] = feature.nCB;
            feature.indCB.push_back(wl);
            if (wl < feature.indCBwl.size()) feature.indCBwl[wl] = feature.nCB;
            feature.nCB++;
        }
    }
    if (feature.nCB==0) return;

    // Compute per-detected-CB counts and total
    std::vector<uint32_t> counts(feature.nCB, 0);
    uint64_t totalRecs = 0;
    for (uint32_t wl=0; wl<perWL.size(); ++wl) {
        uint32_t icb = cbWLtoICB[wl];
        if (icb==(uint32_t)-1) continue;
        counts[icb] = (uint32_t)perWL[wl].size();
        totalRecs += perWL[wl].size();
    }

    // Allocate rGeneUMI and rCBp
    feature.rGeneUMI = new uint32[feature.getRGUStride() * totalRecs];
    feature.rCBp = new uint32*[feature.nCB + 1];
    feature.rCBp[0] = feature.rGeneUMI;
    for (uint32_t i=0; i<feature.nCB; ++i) {
        feature.rCBp[i+1] = feature.rCBp[i] + feature.getRGUStride() * counts[i];
    }

    // Fill triplets per CB in WL order
    for (uint32_t wl=0; wl<perWL.size(); ++wl) {
        uint32_t icb = cbWLtoICB[wl];
        if (icb==(uint32_t)-1) continue;
        uint32_t *blockStart = feature.rCBp[icb];
        uint32_t *dst = blockStart;
        for (const auto &r : perWL[wl]) {
            dst[0] = r.featureId;
            dst[1] = r.umi;
            if (feature.getRGUStride()==3) dst[2] = r.readIndex;
            dst += feature.getRGUStride();
        }
        // Ensure rCBp[icb] points to the start of the block (already set)
        feature.rCBp[icb] = blockStart;
    }

    // Prepare per-CB read counts for caller
    if (feature.rCBn) { delete[] feature.rCBn; feature.rCBn = nullptr; }
    feature.rCBn = new uint32[feature.nCB];
    for (uint32_t iCB=0; iCB<feature.nCB; ++iCB) feature.rCBn[iCB] = counts[iCB];

    // Clear buffers for reuse and free memory
    for (auto &v : perWL) { std::vector<ReadInfoRecord>().swap(v); }
    perWL.clear();
    perWL.shrink_to_fit();
}
