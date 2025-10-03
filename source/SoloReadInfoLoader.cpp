#include "SoloReadInfoLoader.h"
#include "SoloReadFeature.h"
#include "soloInputFeatureUMI.h"
#include "SoloCommon.h"
#include "serviceFuns.cpp"

void SoloReadInfoLoader::load(SoloReadFeature &rf,
                              SoloReadInfoMode mode,
                              const RecordSink &sink,
                              std::vector<uint32_t> &cbReadCountTotal,
                              SoloReadFlagClass &readFlagCounts,
                              std::vector<uint32_t> &nReadPerCBunique1,
                              std::vector<uint32_t> &nReadPerCBmulti1) {

    rf.streamReads->flush();
    rf.streamReads->seekg(0,std::ios::beg);

    uint32 feature;
    uint64 umi, iread, prevIread=(uint64)-1;
    int32 cbmatch;
    int64 cb;
    std::vector<uint32> trIdDist;

    while (soloInputFeatureUMI(rf.streamReads, rf.featureType, rf.readIndexYes, rf.P.sjAll, iread, cbmatch, feature, umi, trIdDist, readFlagCounts)) {
        if (feature == (uint32)(-1) && !rf.readIndexYes) {
            rf.streamReads->ignore((uint32)-1, '\n');
            continue;
        };

        bool readIsCounted = false;
        bool featGood = ( feature != (uint32)(-1) );
        bool noMMtoWLwithoutExact = false;
        bool noTooManyWLmatches = false;

        if (cbmatch<=1) {
            *rf.streamReads >> cb;
            if ( rf.pSolo.CBmatchWL.oneExact && cbmatch==1 && cbReadCountTotal[cb]==0 ) {
                noMMtoWLwithoutExact = true;
                // Missing CB due to oneExact requirement without exact match
                if (sink) {
                    uint32_t umi32 = (uint32_t)-1;
                    uint8_t status = 0;
                    sink({(uint64_t)iread, (uint32_t)-1, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
                }
            } else {
                if (!rf.pSolo.cbWLyes) {
                    cb=binarySearchExact<uintCB>(cb, rf.pSolo.cbWL.data(), rf.pSolo.cbWLsize);
                    if (cb+1 == 0) continue;
                };
                if (featGood) {
                    readIsCounted = true;
                    if (sink) {
                        uint32_t umi32 = (uint32_t)umi;
                        uint8_t status = (umi32==(uint32_t)-1) ? 2 : 1;
                        sink({(uint64_t)iread, (uint32_t)cb, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
                    }
                } else {
                    // no-feature, still emit minimal data; status indicates umi validity
                    if (sink) {
                        uint32_t umi32 = (uint32_t)umi;
                        uint8_t status = (umi32==(uint32_t)-1) ? 2 : 1;
                        sink({(uint64_t)iread, (uint32_t)cb, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
                    }
                };
            };
        } else {
            #ifdef MATCH_CellRanger
            double ptot=0.0, pmax=0.0, pin;
            #else
            float ptot=0.0, pmax=0.0, pin;
            #endif
            for (uint32 ii=0; ii<(uint32)cbmatch; ii++) {
                uint32 cbin; char qin; *rf.streamReads >> cbin >> qin;
                if (cbReadCountTotal[cbin]>0) {
                    qin -= rf.pSolo.QSbase; qin = qin < rf.pSolo.QSmax ? qin : rf.pSolo.QSmax;
                    pin=cbReadCountTotal[cbin]*std::pow(10.0,-qin/10.0);
                    ptot+=pin; if (pin>pmax) { cb=cbin; pmax=pin; };
                };
            };
            if (ptot>0.0 && pmax>=rf.pSolo.cbMinP*ptot) {
                if (featGood) {
                    readIsCounted = true;
                    if (sink) {
                        uint32_t umi32 = (uint32_t)umi;
                        uint8_t status = (umi32==(uint32_t)-1) ? 2 : 1;
                        sink({(uint64_t)iread, (uint32_t)cb, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
                    }
                } else {
                    if (sink) {
                        uint32_t umi32 = (uint32_t)umi;
                        uint8_t status = (umi32==(uint32_t)-1) ? 2 : 1;
                        sink({(uint64_t)iread, (uint32_t)cb, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
                    }
                };
            } else {
                noTooManyWLmatches = true;
                // Emit missing CB status for multi-match failure
                if (sink) {
                    uint32_t umi32 = (uint32_t)-1;
                    uint8_t status = 0;
                    sink({(uint64_t)iread, (uint32_t)-1, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
                }
            };
        };

        if ( !rf.readIndexYes || iread != prevIread ) {
            prevIread = iread;
            if (mode==SoloReadInfoMode::Counting && featGood) {
                if (cbmatch==0) { rf.stats.V[rf.stats.yessubWLmatchExact]++; }
                else if (noMMtoWLwithoutExact) { rf.stats.V[rf.stats.noMMtoWLwithoutExact]++; }
                else if (noTooManyWLmatches) { rf.stats.V[rf.stats.noTooManyWLmatches]++; };
            };
            if (mode==SoloReadInfoMode::Counting && readIsCounted) {
                if (feature<geneMultMark) nReadPerCBunique1[cb]++; else nReadPerCBmulti1[cb]++;
            };
            if ( mode==SoloReadInfoMode::Counting && rf.pSolo.readStatsYes[rf.featureType] ) {
                if ( readIsCounted ) {
                    if ( readFlagCounts.checkBit(readFlagCounts.featureU) ) readFlagCounts.setBit(readFlagCounts.countedU);
                    if ( readFlagCounts.checkBit(readFlagCounts.featureM) ) readFlagCounts.setBit(readFlagCounts.countedM);
                };
                readFlagCounts.setBit(rf.readFlag.cbMatch);
                if (cbmatch==0) { readFlagCounts.setBit(readFlagCounts.cbPerfect); readFlagCounts.countsAdd(cb); }
                else if (cbmatch==1 && !noMMtoWLwithoutExact) { readFlagCounts.setBit(readFlagCounts.cbMMunique); readFlagCounts.countsAdd(cb); }
                else if (cbmatch>1 && !noTooManyWLmatches) { readFlagCounts.setBit(readFlagCounts.cbMMmultiple); readFlagCounts.countsAdd(cb); }
                else { readFlagCounts.countsAddNoCB(); };
            };
        };
    };
}

void SoloReadInfoLoader::loadMinimal(SoloReadFeature &rf, const RecordSink &sink, std::vector<uint32_t> &cbReadCountTotal) {
    SoloReadFlagClass dummyFlags;
    std::vector<uint32_t> nU(cbReadCountTotal.size(), 0), nM(cbReadCountTotal.size(), 0);
    load(rf, SoloReadInfoMode::Minimal, sink, cbReadCountTotal, dummyFlags, nU, nM);
}


