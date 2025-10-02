#include "SoloFeature.h"
#include "streamFuns.h"

SoloFeature::SoloFeature(Parameters &Pin, ReadAlignChunk **RAchunk, Transcriptome &inTrans, int32 feTy, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll)
            : P(Pin), RAchunk(RAchunk), Trans(inTrans), featureType(feTy), soloFeatAll(soloFeatAll), pSolo(P.pSolo), readBarSum(readBarSumIn)
{
    if (featureType>=0) {//otherwise we do not need these arrays - e.g. with --runMode soloCellFiltering 
        readFeatSum = new SoloReadFeature(featureType,P,-1);
        readFeatAll = new SoloReadFeature*[P.runThreadN];
    };
    
    //number of features
    switch (featureType) {
        case SoloFeatureTypes::Gene :
        case SoloFeatureTypes::GeneFull :
        case SoloFeatureTypes::GeneFull_Ex50pAS :
        case SoloFeatureTypes::GeneFull_ExonOverIntron :
        case SoloFeatureTypes::Velocyto :
            featuresNumber=Trans.nGe;
            break;
        case SoloFeatureTypes::SJ :
            featuresNumber=P.sjAll[0].size();
            break;
        default:
            featuresNumber = -1; //undefined
    };    
};

// Shared helper implementations for readInfo management
void SoloFeature::resetPackedStorage(uint32_t nReads)
{
#ifdef SOLO_USE_PACKED_READINFO
    packedReadInfo.init(nReads, pSolo.cbWLstr.size(), pSolo.umiL);
#else
    readInfo.resize(nReads, {(uint64)-1, (uint32)-1});
#endif
}

void SoloFeature::recordReadInfo(uint32_t readId, uint32_t cbIdx, uint32_t umiPacked, uint8_t status)
{
#ifdef SOLO_USE_PACKED_READINFO
    packedReadInfo.set(readId, cbIdx, umiPacked, status);
#else
    readInfo[readId].cb = cbIdx;
    readInfo[readId].umi = umiPacked;
#endif
}

uint32_t SoloFeature::getPackedCB(uint32_t readId) const
{
#ifdef SOLO_USE_PACKED_READINFO
    return packedReadInfo.getCB(readId);
#else
    return (uint32_t)readInfo[readId].cb;
#endif
}

uint32_t SoloFeature::getPackedUMI(uint32_t readId) const
{
#ifdef SOLO_USE_PACKED_READINFO
    return packedReadInfo.getUMI(readId);
#else
    return readInfo[readId].umi;
#endif
}

uint8_t SoloFeature::getPackedStatus(uint32_t readId) const
{
#ifdef SOLO_USE_PACKED_READINFO
    return packedReadInfo.getStatus(readId);
#else
    if (readInfo[readId].cb == (uint64)-1) return 0;
    if (readInfo[readId].umi == (uint32)-1) return 2;
    return 1;
#endif
}

void SoloFeature::clearLarge()
{
    cbFeatureUMImap.clear();
    cbFeatureUMImap.shrink_to_fit();
    countCellGeneUMI.clear();
    countCellGeneUMI.shrink_to_fit();
    countCellGeneUMIindex.clear();
    countCellGeneUMIindex.shrink_to_fit();
    countMatMult.i.clear();
    countMatMult.i.shrink_to_fit();
    countMatMult.m.clear();
    countMatMult.m.shrink_to_fit();
    //indCB.clear(); //needed for Velocyto
    //indCB.shrink_to_fit();
    indCBwl.clear();
    indCBwl.shrink_to_fit();
    nGenePerCB.clear();
    nGenePerCB.shrink_to_fit();
    nGenePerCBmulti.clear();
    nGenePerCBmulti.shrink_to_fit();
    nReadPerCB.clear();
    nReadPerCB.shrink_to_fit();
    nReadPerCBtotal.clear();
    nReadPerCBtotal.shrink_to_fit();
    nReadPerCBunique.clear();
    nReadPerCBunique.shrink_to_fit();
    nUMIperCB.clear();
    nUMIperCB.shrink_to_fit();
    nUMIperCBmulti.clear();
    nUMIperCBmulti.shrink_to_fit();
    nUMIperCBsorted.clear();
    nUMIperCBsorted.shrink_to_fit();
    sjAll[0].clear();
    sjAll[0].shrink_to_fit();
    sjAll[1].clear();
    sjAll[1].shrink_to_fit();
};
