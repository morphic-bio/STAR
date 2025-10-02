#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"

void SoloFeature::prepareReadInfoOnly()
{
    // Minimal processing to populate readInfo without heavy matrix allocations
    // This is called in skipProcessing mode when tag table or CB/UB injection is needed
    
    time_t rawTime;
    
    if (!pSolo.readInfoYes[featureType]) {
        // readInfo not needed for this feature type
        return;
    }
    
    // Allocate and initialize readInfo (same as countCBgeneUMI init)
    resetPackedStorage(nReadsInput);
    time(&rawTime);
#ifdef SOLO_USE_PACKED_READINFO
    P.inOut->logMain << timeMonthDayTime(rawTime) 
                     << " ... Allocated packed readInfo array for skipProcessing mode, nReadsInput = " 
                     << nReadsInput << endl;
#else
    P.inOut->logMain << timeMonthDayTime(rawTime) 
                     << " ... Allocated readInfo array for skipProcessing mode, nReadsInput = " 
                     << nReadsInput << endl;
#endif
    
    // We need to call inputRecords to populate readInfo from the temporary Solo files
    // but we don't need the heavy rGeneUMI array or counting structures
    
    // Minimal setup for inputRecords - we pass dummy arrays where matrices would go
    rguStride = 2;
    if (pSolo.readIndexYes[featureType])
        rguStride = 3;
    
    // Allocate minimal temporary array just for inputRecords to write into
    // This will be thrown away - we only care about readInfo being populated
    rGeneUMI = new uint32[rguStride * nReadsMapped];
    rCBp = new uint32*[nCB + 1];
    uint32 **rCBpa = new uint32*[pSolo.cbWLsize + 1];
    
    rCBp[0] = rGeneUMI;
    rCBpa[0] = rGeneUMI;
    nCB = 0;
    for (uint32 ii = 0; ii < pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii] > 0) {
            rCBp[nCB + 1] = rCBp[nCB] + rguStride * readFeatSum->cbReadCount[ii];
            ++nCB;
        }
        rCBpa[ii + 1] = rCBp[nCB];
    }
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) 
                     << " ... Allocated minimal arrays for skipProcessing " 
                     << nReadsMapped * rguStride * 4.0 / 1024 / 1024 / 1024 << " GiB" << endl;
    
    // Call inputRecords to populate readInfo from temporary files
    readFlagCounts.flagCounts.reserve(nCB * 3 / 2);
    readFlagCounts.flagCountsNoCB = {};
    vector<uint32> nReadPerCBunique1(pSolo.cbWLsize), nReadPerCBmulti1(pSolo.cbWLsize);
    
#ifdef SOLO_USE_PACKED_READINFO
    auto sink = [this](uint64 readId, uint32 cbIdx, uint32 umiPacked, uint8 status){
        recordReadInfo((uint32_t)readId, cbIdx, umiPacked, status);
    };
    for (int ii = 0; ii < P.runThreadN; ii++) {
        readFeatAll[ii]->inputRecords(rCBpa, rguStride, readBarSum->cbReadCountExact, 
                                      readFlagCounts, 
                                      nReadPerCBunique1, nReadPerCBmulti1,
                                      sink);
        readFeatSum->addStats(*readFeatAll[ii]);
    }
#else
    for (int ii = 0; ii < P.runThreadN; ii++) {
        readFeatAll[ii]->inputRecords(rCBpa, rguStride, readBarSum->cbReadCountExact, 
                                      readInfo, readFlagCounts, 
                                      nReadPerCBunique1, nReadPerCBmulti1);
        readFeatSum->addStats(*readFeatAll[ii]);
    }
#endif

    // Instead of manually assigning readInfo here, traverse the rGU structure
    // through collapseUMIall in minimal mode to apply UMI corrections and
    // populate readInfo consistently with the main pipeline.
    countMatStride = 1;
    countCellGeneUMIindex.assign(nCB+1, 0);
    nReadPerCB.resize(nCB);
    nReadPerCBunique.resize(nCB);
    nReadPerCBtotal.resize(nCB);
    nUMIperCB.resize(nCB);
    nGenePerCB.resize(nCB);
    nReadPerCBmax = 0;
    for (uint32 iCB=0; iCB<nCB; iCB++) {
        nReadPerCB[iCB] = (rCBpa[indCB[iCB]] - rCBp[iCB]) / rguStride;
        if (nReadPerCB[iCB] > nReadPerCBmax) nReadPerCBmax = nReadPerCB[iCB];
    }

    collapseUMIall(true);
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) 
                     << " ... Finished populating readInfo from Solo files (skipProcessing mode)" << endl;
    
    // Clean up the temporary arrays immediately - we don't need them
    delete[] rGeneUMI;
    delete[] rCBp;
    delete[] rCBpa;
    
    rGeneUMI = nullptr;
    rCBp = nullptr;
}

