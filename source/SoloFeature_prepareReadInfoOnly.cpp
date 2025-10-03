#include "SoloFeature.h"
#include "SoloReadInfoLoader.h"
#include "SoloReadInfoSink.h"
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
    P.inOut->logMain << timeMonthDayTime(rawTime) 
                     << " ... Allocated packed readInfo array for skipProcessing mode, nReadsInput = " 
                     << nReadsInput << endl;
    
    // Use loader + MinimalSink to populate readInfo directly without allocating rGeneUMI/CB arrays
    readFlagCounts.flagCounts.clear();
    readFlagCounts.flagCountsNoCB = {};

    SoloReadInfoLoader loader;
    MinimalSink sink;
    for (int ii = 0; ii < P.runThreadN; ii++) {
        loader.loadMinimal(*readFeatAll[ii],
                           [&](const ReadInfoRecord &rec){ sink.onRecord(*this, rec); },
                           readBarSum->cbReadCountExact);
        readFeatSum->addStats(*readFeatAll[ii]);
    }
    sink.finalize(*this);

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Finished populating readInfo via loader/minimal sink (skipProcessing mode)" << endl;
}

