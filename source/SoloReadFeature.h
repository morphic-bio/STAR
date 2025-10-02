#ifndef H_SoloReadFeature
#define H_SoloReadFeature
#include <set>
#include <map>
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "SoloReadBarcode.h"
#include "SoloCommon.h"
#include "SoloReadFeatureStats.h"
#include "ReadAnnotations.h"
#include <functional>

class SoloFeature;

class SoloReadFeature {
public:

    uint32 homoPolymer[4];//homopolymer constants

    vector<uint32> cbReadCount;
    map <uintCB,uint32> cbReadCountMap;
    
    vector<uint32> transcriptDistCount;
    
    bool readInfoYes ,readIndexYes;

    fstream *streamReads;

    string cbSeq, umiSeq, cbQual, umiQual;

    SoloReadFlagClass readFlag;
    
    SoloReadFeatureStats stats;

    SoloReadFeature (int32 feTy, Parameters &Pin, int iChunk);
    void record(SoloReadBarcode &soloBar, uint nTr, Transcript **alignOut, uint64 iRead, ReadAnnotations &readAnnot);
    void addCounts(const SoloReadFeature &soloCBin);
    void addStats(const SoloReadFeature &soloCBin);
    void statsOut(ofstream &streamOut);
    void inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, vector<readInfoStruct> &readInfo, SoloReadFlagClass &readFlagCounts,
                      vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1);
    // Overload that emits read info via a sink (avoids requiring legacy vector storage)
    void inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, SoloReadFlagClass &readFlagCounts,
                      vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1,
                      const std::function<void(uint64, uint32, uint32, uint8)> &recordSink);

private:
    const int32 featureType;

    Parameters &P;
    ParametersSolo &pSolo;
};

#endif
