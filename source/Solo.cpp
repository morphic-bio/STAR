#include "Solo.h"
#include "TimeFunctions.h"
#include "streamFuns.h"
#include "BAMTagBuffer.h"

Solo::Solo(ReadAlignChunk **RAchunkIn, Parameters &Pin, Transcriptome &inTrans)
                       :  RAchunk(RAchunkIn), P(Pin), Trans(inTrans), pSolo(P.pSolo)
{
    if ( pSolo.type == 0 )
        return;
    
    readBarSum = new SoloReadBarcode(P);
    
    if ( pSolo.type == pSolo.SoloTypes::CB_samTagOut )
        return;

    // Note: soloFeat allocation happens regardless of skipProcessing
    // The skip flag only affects post-mapping processing in processAndOutput()
    soloFeat = new SoloFeature*[pSolo.nFeatures];
    for (uint32 ii=0; ii<pSolo.nFeatures; ii++)
        soloFeat[ii] = new SoloFeature(P, RAchunk, Trans, pSolo.features[ii], readBarSum, soloFeat);
};

///////////////////////////////////////////////////////////////////////////////////// post-mapping processing only
//overloaded: only soloCellFiltering
Solo::Solo(Parameters &Pin, Transcriptome &inTrans)
          :  P(Pin), Trans(inTrans), pSolo(P.pSolo)
{
    if ( P.runMode != "soloCellFiltering" )
        return; //passing through, return back to executing STAR
        
    time_t timeCurrent;

    time( &timeCurrent);
    *P.inOut->logStdOut << timeMonthDayTime(timeCurrent) << " ..... starting SoloCellFiltering" <<endl;
    
    soloFeat = new SoloFeature*[1];
    
    soloFeat[0] = new SoloFeature(P, NULL, Trans, -1, NULL, soloFeat);
    soloFeat[0]->loadRawMatrix();
    soloFeat[0]->cellFiltering();
    
    time( &timeCurrent);
    *P.inOut->logStdOut << timeMonthDayTime(timeCurrent) << " ..... finished successfully\n" <<flush;
    P.inOut->logMain  << "ALL DONE!\n" << flush;
    exit(0);
};


////////////////////////////////////////////////////////////////////////////////////
void Solo::processAndOutput()
{
    if (pSolo.type==0 )
        return;
    
    {//collect barcode statistics    
        if (pSolo.cbWLyes) {//now we can define WL and counts
            for (int ii=0; ii<P.runThreadN; ii++) {
                readBarSum->addCounts(*RAchunk[ii]->RA->soloRead->readBar);
                readBarSum->addStats(*RAchunk[ii]->RA->soloRead->readBar);
                delete RAchunk[ii]->RA->soloRead->readBar; //not needed anymore
            };
        };

        ofstream *statsStream = &ofstrOpen(P.outFileNamePrefix+pSolo.outFileNames[0]+"Barcodes.stats",ERROR_OUT, P);
        readBarSum->statsOut(*statsStream);
        statsStream->close();

        //pseudocounts
        if (pSolo.CBmatchWL.mm1_multi_pc) {
            for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
                readBarSum->cbReadCountExact[ii]++;//add one to exact counts
            };
        };
    };
    
    if (pSolo.type==pSolo.SoloTypes::CB_samTagOut)
        return;

    // Early exit for skipProcessing mode with minimal readInfo preparation
    if (pSolo.skipProcessing) {
        // If tag table or CB/UB injection is enabled, we need minimal readInfo population
        if (pSolo.writeTagTableEnabled || pSolo.addTagsToUnsorted) {
            // Process only the samAttrFeature to populate readInfo without full counting
            uint32 attrFeatureIdx = pSolo.featureInd[pSolo.samAttrFeature];
            if (attrFeatureIdx < pSolo.nFeatures) {
                *P.inOut->logStdOut << timeMonthDayTime() << " ..... preparing readInfo for tag export (skipProcessing mode)\n" <<flush;
                P.inOut->logMain    << timeMonthDayTime() << " ..... preparing readInfo for tag export (skipProcessing mode)\n" <<flush;
                
                // Call processRecords which will populate readInfo
                // The feature processing code should handle skipProcessing internally
                soloFeat[attrFeatureIdx]->processRecords();
            }
        }
        
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... skipping Solo counting (soloSkipProcessing=yes)\n" <<flush;
        P.inOut->logMain    << timeMonthDayTime() << " ..... skipping Solo counting (soloSkipProcessing=yes)\n" <<flush;
        return;
    }

    {//process all features
        *P.inOut->logStdOut << timeMonthDayTime() << " ..... started Solo counting\n" <<flush;
        P.inOut->logMain    << timeMonthDayTime() << " ..... started Solo counting\n" <<flush;

        for (uint32 ii=0; ii<pSolo.nFeatures; ii++) {
            soloFeat[ii]->processRecords();
            //if (!pSolo.readInfoYes[soloFeat[ii]->featureType]) {
            //    delete soloFeat[ii];
            //};
        };

        *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished Solo counting\n" <<flush;
        P.inOut->logMain    << timeMonthDayTime() << " ..... finished Solo counting\n" <<flush;
    };
};
